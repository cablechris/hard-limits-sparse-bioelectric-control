from __future__ import annotations


def edge_key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def lattice_index(row: int, col: int, size: int) -> int:
    return row * size + col


def build_square_neighbors(size: int) -> dict[int, tuple[int, ...]]:
    neighbors: dict[int, tuple[int, ...]] = {}
    for row in range(size):
        for col in range(size):
            site = lattice_index(row, col, size)
            local: list[int] = []
            if row > 0:
                local.append(lattice_index(row - 1, col, size))
            if row + 1 < size:
                local.append(lattice_index(row + 1, col, size))
            if col > 0:
                local.append(lattice_index(row, col - 1, size))
            if col + 1 < size:
                local.append(lattice_index(row, col + 1, size))
            neighbors[site] = tuple(local)
    return neighbors
