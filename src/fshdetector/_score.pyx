# distutils: language = c
"""Compares and combines unit-barcode matrices."""

from libc.stdint cimport int8_t, int32_t

import numpy as np

cimport numpy as cnp

cnp.import_array()


def score_offsets(
    const int8_t[:, ::1] template,
    const int8_t[:, ::1] query,
    const double[:, :, ::1] table,
    const double[::1] noise,
    const long[::1] offsets,
):
    """Compares a query with a template at each offset.

    Row ``j`` of the query is paired with row ``j + offset`` of the template.
    A site counts only if both units carry an observation (value >= 0).

    Args:
        template: units x sites, allele codes, -1 where unobserved.
        query: units x sites.
        table: site x allele x allele log likelihood ratios.
        noise: per site probability that two reads of one unit disagree.
        offsets: placements to test.

    Returns:
        A tuple of arrays, one value per offset: score, compared sites,
        mismatched sites, expected mismatches, variance of that expectation,
        and units that share at least one site.
    """
    cdef Py_ssize_t n_offsets = offsets.shape[0]
    cdef Py_ssize_t n_template = template.shape[0]
    cdef Py_ssize_t n_query = query.shape[0]
    cdef Py_ssize_t n_sites = template.shape[1]
    if query.shape[1] != n_sites or table.shape[0] != n_sites:
        raise ValueError("template, query, and table disagree on sites")

    score = np.zeros(n_offsets)
    sites = np.zeros(n_offsets, dtype=np.int64)
    mismatches = np.zeros(n_offsets, dtype=np.int64)
    expected = np.zeros(n_offsets)
    variance = np.zeros(n_offsets)
    units = np.zeros(n_offsets, dtype=np.int64)
    cdef double[::1] v_score = score
    cdef long[::1] v_sites = sites
    cdef long[::1] v_mismatches = mismatches
    cdef double[::1] v_expected = expected
    cdef double[::1] v_variance = variance
    cdef long[::1] v_units = units

    cdef Py_ssize_t i, j, row, site
    cdef long offset
    cdef int8_t a, b
    cdef bint shared
    with nogil:
        for i in range(n_offsets):
            offset = offsets[i]
            for j in range(n_query):
                row = j + offset
                if row < 0 or row >= n_template:
                    continue
                shared = False
                for site in range(n_sites):
                    a = template[row, site]
                    b = query[j, site]
                    if a < 0 or b < 0:
                        continue
                    shared = True
                    v_score[i] += table[site, a, b]
                    v_sites[i] += 1
                    v_mismatches[i] += a != b
                    v_expected[i] += noise[site]
                    v_variance[i] += noise[site] * (1.0 - noise[site])
                v_units[i] += shared
    return score, sites, mismatches, expected, variance, units


def score_units(
    const int8_t[:, ::1] template,
    const int8_t[:, ::1] query,
    const double[:, :, ::1] table,
    const double[::1] noise,
):
    """Compares every query unit with every template unit.

    Args:
        template: units x sites, allele codes, -1 where unobserved.
        query: units x sites.
        table: site x allele x allele log likelihood ratios.
        noise: per site probability that two reads of one unit disagree.

    Returns:
        A tuple of template-units x query-units arrays: score, compared sites,
        mismatched sites, expected mismatches, and variance of that
        expectation.
    """
    cdef Py_ssize_t n_template = template.shape[0]
    cdef Py_ssize_t n_query = query.shape[0]
    cdef Py_ssize_t n_sites = template.shape[1]
    if query.shape[1] != n_sites or table.shape[0] != n_sites:
        raise ValueError("template, query, and table disagree on sites")

    shape = (n_template, n_query)
    score = np.zeros(shape)
    sites = np.zeros(shape, dtype=np.int64)
    mismatches = np.zeros(shape, dtype=np.int64)
    expected = np.zeros(shape)
    variance = np.zeros(shape)
    cdef double[:, ::1] v_score = score
    cdef long[:, ::1] v_sites = sites
    cdef long[:, ::1] v_mismatches = mismatches
    cdef double[:, ::1] v_expected = expected
    cdef double[:, ::1] v_variance = variance

    cdef Py_ssize_t i, j, site
    cdef int8_t a, b
    with nogil:
        for i in range(n_template):
            for j in range(n_query):
                for site in range(n_sites):
                    a = template[i, site]
                    b = query[j, site]
                    if a < 0 or b < 0:
                        continue
                    v_score[i, j] += table[site, a, b]
                    v_sites[i, j] += 1
                    v_mismatches[i, j] += a != b
                    v_expected[i, j] += noise[site]
                    v_variance[i, j] += noise[site] * (1.0 - noise[site])
    return score, sites, mismatches, expected, variance


def add_counts(
    int32_t[:, :, ::1] counts, const int8_t[:, ::1] matrix, Py_ssize_t offset
):
    """Adds the alleles of a read to a units x sites x alleles count table.

    Row ``j`` of the read goes to row ``j + offset`` of the table; rows outside
    the table are ignored.
    """
    if matrix.shape[1] != counts.shape[1]:
        raise ValueError("matrix and counts disagree on sites")
    cdef Py_ssize_t j, row, site
    cdef int8_t a
    with nogil:
        for j in range(matrix.shape[0]):
            row = j + offset
            if row < 0 or row >= counts.shape[0]:
                continue
            for site in range(matrix.shape[1]):
                a = matrix[j, site]
                if 0 <= a < counts.shape[2]:
                    counts[row, site, a] += 1


def consensus(const int32_t[:, :, ::1] counts, Py_ssize_t rows):
    """The most common allele per unit and site; ties and no data are -1."""
    result = np.full((rows, counts.shape[1]), -1, dtype=np.int8)
    cdef int8_t[:, ::1] v_result = result
    cdef Py_ssize_t row, site, allele
    cdef int32_t best
    cdef int8_t winner
    with nogil:
        for row in range(min(rows, counts.shape[0])):
            for site in range(counts.shape[1]):
                best = 0
                winner = -1
                for allele in range(counts.shape[2]):
                    if counts[row, site, allele] > best:
                        best = counts[row, site, allele]
                        winner = <int8_t>allele
                    elif counts[row, site, allele] == best:
                        winner = -1
                v_result[row, site] = winner
    return result
