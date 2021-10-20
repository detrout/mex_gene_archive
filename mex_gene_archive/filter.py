#######
# Functions for making filtered .mtx file
import csv


def read_barcode_lineno_map(stream):
    """Build a map of barcodes to line number from a stream"""
    barcodes = {}
    reader = csv.reader(stream, delimiter="\t")
    for i, line in enumerate(reader):
        barcodes[line[0]] = i + 1

    return barcodes


def compute_raw_to_filtered_map(filtered_barcodes, raw_barcodes):
    """Generate a raw index to filtered index mapping"""
    raw_to_filtered_mapping = {}
    for filtered_barcode in filtered_barcodes:
        filtered_index = filtered_barcodes[filtered_barcode]
        raw_index = raw_barcodes[filtered_barcode]
        raw_to_filtered_mapping[raw_index] = filtered_index
    return raw_to_filtered_mapping


def filter_mtx(raw_barcode_filename, raw_matrix_filename, filtered_barcode_filename):
    """generating a filtered mtx file from raw and filtered values

    We read the the raw barcodes, raw matrix, and filtered barcodes
    and then yield the rows from the raw matrix file to only include
    the barcodes in the filtered barcode list and rewrite the index values
    to match the filtered barcode list.
    """
    with open(raw_barcode_filename, "rt") as instream:
        raw_barcodes = read_barcode_lineno_map(instream)
    with open(filtered_barcode_filename, "rt") as instream:
        filtered_barcodes = read_barcode_lineno_map(instream)
    raw_to_filtered_mapping = compute_raw_to_filtered_map(
        filtered_barcodes, raw_barcodes
    )

    header = True
    results = []
    with open(raw_matrix_filename, "rt") as instream:
        # copy comments
        for line in instream:
            if line.startswith("%"):
                yield line
            elif header:
                # After the comment comes the one header line
                total_rows, total_columns, total_counts = [
                    int(x) for x in line.rstrip().split()
                ]
                assert total_columns == len(raw_barcodes)
                header = False
            else:
                # row, column, count
                row, column, count = line.rstrip().split()
                row = int(row)
                column = int(column)
                if column in raw_to_filtered_mapping:
                    new_column = raw_to_filtered_mapping[column]
                    results.append((row, new_column, count))

        rs = sorted(results, key=lambda row: (row[1], row[0]))
        total_columns = len(filtered_barcodes)
        total_counts = len(rs)
        yield "{} {} {}\n".format(total_rows, total_columns, total_counts)
        for row, column, count in rs:
            yield "{} {} {}\n".format(row, column, count)


def write_filtered_mtx(
    raw_barcode_filename,
    raw_matrix_filename,
    filtered_barcode_filename,
    filtered_matrix_filename,
):
    with open(filtered_matrix_filename, "wt") as outstream:
        for line in filter_mtx(
            raw_barcode_filename, raw_matrix_filename, filtered_barcode_filename
        ):
            outstream.write(line)
