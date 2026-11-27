# distutils: language = c
"""Exact k-mer screen that selects reads before alignment.

A read is kept if it shares at least one k-mer with a query sequence on either
strand. Extra sequences can add k-mers; their low-complexity k-mers are left
out, because simple repeats such as poly-A occur in many unrelated reads. Only every `stride`-th read position is tested. The screen cannot lose
a read that has an exact match of `k + stride - 1` bases or more, because one
of the tested k-mers then lies fully inside that match.
"""

from libc.stdint cimport uint64_t, uint8_t
from libc.stdlib cimport free, malloc

cdef uint64_t EMPTY = <uint64_t>0xFFFFFFFFFFFFFFFF

cdef uint8_t CODE[256]


cdef void _fill_codes():
    cdef int i
    cdef const unsigned char *upper = b"ACGT"
    cdef const unsigned char *lower = b"acgt"
    for i in range(256):
        CODE[i] = 4
    for i in range(4):
        CODE[upper[i]] = i
        CODE[lower[i]] = i


_fill_codes()


cdef int _MIN_DINUCLEOTIDES = 6


cdef bint _complex(uint64_t key, int k) noexcept:
    """True if the k-mer holds enough different dinucleotides."""
    cdef int i, found = 0
    cdef uint64_t seen = 0, bit
    for i in range(k - 1):
        bit = <uint64_t>1 << ((key >> (2 * i)) & 15)
        if not seen & bit:
            seen |= bit
            found += 1
    return found >= _MIN_DINUCLEOTIDES


cdef inline uint64_t _hash(uint64_t key) noexcept nogil:
    key ^= key >> 31
    key *= <uint64_t>0x9E3779B97F4A7C15
    return key ^ (key >> 29)


cdef class KmerScreen:
    """Hash set of query k-mers.

    Args:
        sequences: query sequences; both strands are indexed.
        k: k-mer length, 1 to 31.
        stride: distance between tested read positions.
        extra: more sequences, of which only complex k-mers are indexed.
    """

    cdef uint64_t *table
    cdef uint64_t mask
    cdef readonly int k
    cdef readonly int stride

    def __cinit__(self, sequences, int k=21, int stride=8, extra=()):
        if not 1 <= k <= 31:
            raise ValueError("k must be between 1 and 31")
        if stride < 1:
            raise ValueError("stride must be at least 1")
        self.k = k
        self.stride = stride
        keys = set()
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        for sequence in sequences:
            for strand in (sequence, sequence.translate(complement)[::-1]):
                keys.update(_kmers(strand.encode("ascii"), k))
        for sequence in extra:
            for strand in (sequence, sequence.translate(complement)[::-1]):
                keys.update(
                    key
                    for key in _kmers(strand.encode("ascii"), k)
                    if _complex(key, k)
                )
        cdef uint64_t size = 1024
        while size < 4 * <uint64_t>len(keys):
            size *= 2
        self.table = <uint64_t *>malloc(size * sizeof(uint64_t))
        if self.table is NULL:
            raise MemoryError()
        self.mask = size - 1
        cdef uint64_t i, key, slot
        for i in range(size):
            self.table[i] = EMPTY
        for key in keys:
            slot = _hash(key) & self.mask
            while self.table[slot] != EMPTY:
                slot = (slot + 1) & self.mask
            self.table[slot] = key

    def __dealloc__(self):
        free(self.table)

    cpdef bint matches(self, str sequence):
        """Returns True if a tested read k-mer is in the query set."""
        cdef bytes raw = sequence.encode("ascii")
        cdef const unsigned char *data = raw
        cdef Py_ssize_t length = len(raw)
        cdef Py_ssize_t position = 0
        cdef Py_ssize_t offset
        cdef uint64_t key, slot
        cdef uint8_t code
        cdef bint valid
        cdef bint found = False
        with nogil:
            while position + self.k <= length and not found:
                key = 0
                valid = True
                for offset in range(self.k):
                    code = CODE[data[position + offset]]
                    if code > 3:
                        valid = False
                        break
                    key = (key << 2) | code
                if valid:
                    slot = _hash(key) & self.mask
                    while self.table[slot] != EMPTY:
                        if self.table[slot] == key:
                            found = True
                            break
                        slot = (slot + 1) & self.mask
                position += self.stride
        return found


def _kmers(bytes sequence, int k):
    """Yields the 2-bit code of every k-mer without ambiguous bases."""
    cdef uint64_t key = 0
    cdef uint64_t keep = (<uint64_t>1 << (2 * k)) - 1
    cdef int run = 0
    cdef uint8_t code
    cdef const unsigned char *data = sequence
    cdef Py_ssize_t i
    for i in range(len(sequence)):
        code = CODE[data[i]]
        if code > 3:
            run = 0
            key = 0
            continue
        key = ((key << 2) | code) & keep
        run += 1
        if run >= k:
            yield key
