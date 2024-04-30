"""
Program for consolidating a csv file of sequence counts based on miRNA name.
To adjust which files the program attempts to consol, change what glob.glob
searches for at the top of the main() function.
To adjust the consol output file names, change the 2nd arg in filename.replace()
on the first line of the write_file() function.
"""
import csv
import re
import glob


def read_file(filename: str) -> tuple[list, list[list]]:
    '''
    Reads miRNA count data from the given csv file and generates a nested list.

    Args:
        filename: the name of the miRNA count csv file to be consolidated
    Returns:
        header: a list of column lables (str) from the miRNA raw count file
        data: a 2D list of miRNA names (str) and counts (int)
    '''
    with open(filename, 'r', newline='', encoding='utf-8') as in_file:
        reader = csv.reader(in_file)
        data = []
        for line in reader:
            data.append(line)
        header = data.pop(0)
        for row in data:
            # written to skip 1st column (str miRNA names)
            for indx, count in enumerate(row[1:]):
                row[indx+1] = int(count)
        return header, data


def generate_consol_names(data: list[list]) -> list[list]:
    '''
    Appends a consolidatable miRNA name to the end of each data row line.

    Args:
        data: a 2D list of miRNA names (str) and counts (int)
    Returns:
        data: a 2D list of full miRNA names (str), counts (int), and consol names (str)
    '''
    for row in data:
        name_pieces = re.split(r'\W+', row[0])          # splits name on special characters
        if name_pieces[1].lower() == 'mir':             # handles majority of names
            consol_name = name_pieces[2]
            if len(name_pieces) == 5:                   # adds subtype if applicable (ex: 281-2)
                consol_name += '-' + name_pieces[3]
        else:
            consol_name = name_pieces[1]                # handles non-numeric names (ex: bantam)
            if consol_name[:3].lower() == 'mir':
                consol_name = consol_name[3:]
            elif consol_name.lower() == 'let':
                consol_name += '-' + name_pieces[2]
        row.append(consol_name)
    return data


def sort_and_sum(data:list[list]) -> list[list]:
    '''
    Reorders the data rows by their consolidated names, 
    which moves rows with identical consol names adjacent to each other.
    Compiles the counts of rows with identical consol names.

    Args:
        data: a 2D list of full miRNA names (str), counts (int), and consol names (str)
    Returns:
        data: a sorted, summed 2D list of consol miRNA names (str) and counts (int)
    '''
    # sort by name
    data.sort(key = lambda x: x[-1])
    #sum by matching names
    for indx, line in enumerate(data):
        if indx != len(data)-1:
            next_line = data[indx+1]

            if line[-1] == next_line[-1]:
                for i in range(len(line)-2):
                    # actually sum the counts
                    next_line[i+1] += line[i+1]
    # delete already summed rows
    indxx = 0
    while indxx < len(data)-1:
        if data[indxx][-1] == data[indxx+1][-1]:
            data.pop(indxx)
        else:
            indxx += 1
    # clean up names
    for line in data:
        line[0] = line.pop()
    return data


def write_file(filename: str, header: list, data: list[list]):
    '''
    Writes the contents of the header & consolidated data to a csv file.

    Args:
        filename: the name of the miRNA raw count file, the basis for the consol file name
        header: a list of column lables (str) from the miRNA raw count file
        data: a sorted, summed 2D list of consol miRNA names (str) and counts (int)
    '''
    output_name = filename.replace('raw_copy', 'py_consol')
    with open(output_name, 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(header)
        for line in data:
            writer.writerow(line)


def main():
    '''
    Consolidates the miRNA count data from each file glob.glob finds,
    and saves the result as a csv file based on the name of the glob input.
    If a consol csv file with that name already exists, it is overwritten.
    '''
    for file_name in glob.glob('**/*raw_copy.csv', recursive=True):
        column_headers, mirna_data = read_file(file_name)
        consol_names = generate_consol_names(mirna_data)
        sorted_data = sort_and_sum(consol_names)
        write_file(file_name, column_headers, sorted_data)


if __name__ == '__main__':
    main()
