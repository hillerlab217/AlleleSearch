import os
import csv
import argparse

"""
adjusted alignlen to the actrrual length without gap
adjusted the subject end to actrual end
made the script runable with one line
"""

# Run example
# python code/blastParser.py --directory output/blast --querys 1 --blast_type blastn --max_Eval 10 --sequence_identity 0 --min_length 0 --file_name output/blast_summary_all.csv

def is_tabular_fmt(lines):
    # Heuristic: a non-empty, non-comment line with >=12 tab-separated fields and no "Query=" token
    for L in lines:
        if not L.strip() or L.startswith('#'):
            continue
        if '\t' in L and len(L.rstrip('\n').split('\t')) >= 12 and "Query=" not in L:
            return True
        break
    return False

def parseDir(directory, querys, blast_type, max_Eval, sequence_identity, min_length):
    listOfInputs = []
    for file in os.listdir(directory):
        if file.endswith(".out"):
            listOfInputs.append(os.path.join(directory, file))

    dataParsed = dict()
    for file in listOfInputs:
        with open(file, "r") as inputFile:
            inputName = os.path.basename(file).split('.')[0]
            inputString = inputFile.readlines()

        if querys == 1:
            if is_tabular_fmt(inputString):
                # implement parseBLAST_single_query_tabular similarly, grouping by contig
                dataParsed = parseBLAST_single_query_tabular(inputName, inputString, dataParsed,
                                                             blast_type, max_Eval, sequence_identity, min_length)
            else:
                dataParsed = parseBLAST(inputName, inputString, dataParsed,
                                        blast_type, max_Eval, sequence_identity, min_length)
        else:
            if is_tabular_fmt(inputString):
                dataParsed = parseBLAST_multi_query_tabular(inputName, inputString, dataParsed,
                                                            blast_type, max_Eval, sequence_identity, min_length)
            else:
                dataParsed = parseBLAST_multi_query(inputName, inputString, dataParsed,
                                                    blast_type, max_Eval, sequence_identity, min_length)
    return dataParsed


#parseBLAST takes the contents of the BLAST file and outputs a dictionary
#which contains ID, contig, eval, alignment length and sequence identity, query start, subject start and end
def parseBLAST(entryName, inputString, dataParsed, blast_type, max_Eval, sequence_identity, min_length):
    out = dict()
    contigParsedLines = [[]] #create a list of lists currently size 1
    i = 0
    for line in inputString: #separate by contig
        if ">" in line:    
            i += 1
            contigParsedLines.append([]) #append an empty list
        contigParsedLines[i].append(line) #contains list of from one contig to the next
            
    match = [[]]
    b = 0

    #separate by matches
    for entry in contigParsedLines[1:]:
        if entry[0].startswith('>'):
            contig = entry[0][1:].strip().split()[0]  # gets the contig ID after '>'
        else:
            contig = "UNKNOWN"

        for line in entry:
            if "Score" in line:
                b = b + 1
                match.append([])
                match[b].append(contig) #subject string name
            match[b].append(line)

    #for each match
    for item in match[1:]:
        contig = item[0]
        if blast_type == 'blastn':
            Eval = float(item[1].split("=")[-1].strip())
            strand = item[3].split("/")[-1] #reading frame
        elif blast_type == 'tblastn':
            Eval = float(item[1].split("=")[-1].split()[0].strip(','))
            strand = item[3].split("=")[-1]
        identity = int(item[2].split("/")[1].split("(")[1].split("%")[0]) #sequence identity
        queryStart = int(item[5].split()[1]) #start position of query
        subStart = int(item[7].split()[1]) #start position of subject
        sequences = item[7:]
        seq = ''
        for line in sequences:
            if 'Sbjct' in line:
                seq += line.split()[2].replace('-', '')
        alignlen = len(seq)

        if blast_type == 'blastn':
            if "Minus" in strand:
                subEnd = subStart - alignlen + 1 #keep in nucleotide sequence
            else:
                subEnd = subStart + alignlen - 1
        elif blast_type == 'tblastn':
            if "-" in strand:
                subEnd = subStart - (alignlen * 3) + 1 #keep in nucleotide sequence
            else:
                subEnd = subStart + (alignlen * 3) - 1

        if Eval < max_Eval and identity >= sequence_identity and alignlen >= min_length:
            if contig in out:
                out[contig].append((Eval, alignlen, identity, queryStart, subStart, subEnd, seq))
            else:
                out[contig] = []
                out[contig].append((Eval, alignlen, identity, queryStart, subStart, subEnd, seq))

    dataParsed[entryName] = out
    return dataParsed

def parseBLAST_single_query_tabular(entryName, inputString, dataParsed,
                                    blast_type, max_Eval, sequence_identity, min_length):
    """
    Parse BLAST outfmt 6 (tabular) for the single-query case.

    Columns (first 12 expected):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    Produces:
      dataParsed[entryName][<query_key>] = list of
        (contig, Eval, alignlen, identity, queryStart, subStart, subEnd, seq)
    """
    hits_by_q = {}

    for raw in inputString:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue

        cols = line.split('\t')
        if len(cols) < 12:
            continue

        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = cols[:12]

        # Cast numeric fields
        try:
            Eval       = float(evalue.replace('e', 'E'))
            alignlen   = int(length)
            identity   = float(pident)
            qstart     = int(qstart)
            sstart     = int(sstart)
            send       = int(send)
        except ValueError:
            continue

        # Apply thresholds (min_length treated as absolute bases/aa)
        if (Eval < max_Eval) and (identity >= sequence_identity) and (alignlen >= min_length):
            seq = ""  # outfmt 6 has no aligned sequence text
            hits_by_q.setdefault(qseqid, []).append(
                (sseqid, Eval, alignlen, identity, qstart, sstart, send, seq)
            )

    # For single-query workflows, many downstream writers expect the inner key
    # to be the "entryName". If we only saw one qseqid, map it to entryName.
    if len(hits_by_q) == 1:
        only_q = next(iter(hits_by_q))
        dataParsed[entryName] = {entryName: hits_by_q[only_q]}
    else:
        # Fallback: keep per-qseqid grouping (still valid for most writers)
        dataParsed[entryName] = hits_by_q

    return dataParsed


def parseBLAST_multi_query_tabular(entryName, inputString, dataParsed, blast_type, max_Eval, sequence_identity, min_length):
    """
    Parse BLAST outfmt 6 (tabular):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    Produces: dataParsed[entryName][query] = list of
      (contig, Eval, alignlen, identity, queryStart, subStart, subEnd, seq)
    """
    out = {}
    for raw in inputString:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if len(cols) < 12:
            continue

        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = cols[:12]
        # numeric casts
        try:
            Eval = float(evalue.replace('e','E'))  # handle sci notation
        except ValueError:
            continue
        alignlen   = int(length)
        identity   = float(pident)
        qstart     = int(qstart)
        sstart     = int(sstart)
        send       = int(send)

        # thresholds (treat min_length as an absolute length here)
        if Eval < max_Eval and identity >= sequence_identity and alignlen >= min_length:
            seq = ""  # outfmt 6 doesn't include aligned sequence text
            out.setdefault(qseqid, []).append(
                (sseqid, Eval, alignlen, identity, qstart, sstart, send, seq)
            )

    dataParsed[entryName] = out
    return dataParsed


#same as parseBlast but with multiple queries
#final dictionary contains:
# ID, query, contig, eval, alignment length and sequence identity, query start, subject start and end
def parseBLAST_multi_query(entryName, inputString, dataParsed, blast_type, max_Eval, sequence_identity, min_length):
    out = dict()
    queryParsedLines = [[]] #create a list of lists currently size 1
    i = 0
    for line in inputString: #separate the queries
        if "Query=" in line:
            i += 1
            queryParsedLines.append([]) #append an empty list
        queryParsedLines[i].append(line) #contains list of from one Query= to the next Query=

    entryParsedLines = dict() #[q:entries,qlen]
    for queryGroup in queryParsedLines[1:]: #for each Query=
        entries = [[]]
        v = 0
        q = (queryGroup[0].split("= ")[1] + " " + queryGroup[1].strip()).strip("Query= ").strip()
        #separate by contig for each query
        for line in queryGroup[3:]:
            if ">" in line:
                v = v + 1
                entries.append([])
            entries[v].append(line) #contains information of only one contig
        entryParsedLines[q] = (entries)

    #separate by matches
    match = [[]]
    b=0
    for query in entryParsedLines.keys(): #for each query
        entries = entryParsedLines[query]
        for entry in entries[1:]:
            if entry != []:
                label = entry[0].split()[1:]
                contig = ' '
                for i in label:
                    contig += i + ' '
                contig = contig.strip()
            for line in entry:
                if "Score" in line:
                    b = b + 1
                    match.append([])
                    match[b].append(contig) #subject string name
                    match[b].append(query) #query peptide
                match[b].append(line)

    for item in match[1:]:
        contig = item[0]
        query = item[1]
        if blast_type == 'blastn':
            Eval = float(item[2].split("=")[-1].strip())
            strand = item[4].split("/")[-1] #reading frame
        elif blast_type == 'tblastn':
            Eval = float(item[2].split("=")[-1].split()[0].strip(','))
            strand = item[4].split("=")[-1]
        identity = int(item[3].split("/")[1].split("(")[1].split("%")[0]) #sequence identity
        subStart = int(item[8].split()[1]) #starting position of subject
        sequences = item[8:]
        seq = ''
        for line in sequences:
            if 'Sbjct' in line:
                seq += line.split()[2].replace('-', '')
        alignlen = len(seq)
        queryStart = int(item[6].split()[1]) #starting position of query

        if blast_type == 'blastn':
            if "Minus" in strand:
                subEnd = subStart - alignlen + 1 #keep in nucleotide sequence
            else:
                subEnd = subStart + alignlen - 1
        elif blast_type == 'tblastn':
            if "-" in strand:
                subEnd = subStart - (alignlen * 3) + 1 #keep in nucleotide sequence
            else:
                subEnd = subStart + (alignlen * 3) - 1

        if Eval < max_Eval and identity >= sequence_identity and alignlen >= min_length:
            if query in out:
                out[query].append((contig, Eval, alignlen, identity, queryStart, subStart, subEnd, seq))
            else:
                out[query] = []
                out[query].append((contig, Eval, alignlen, identity, queryStart, subStart, subEnd, seq))

    dataParsed[entryName] = out
    return dataParsed

def output_csv(output, querys, file_name):
    file_name = file_name.strip() + '.csv'
    
    with open(file_name, mode='w', newline='') as file:
        writer = csv.writer(file)

        if querys == 1:
            # Write header
            titles = ['ID', 'Contig', 'E-value', 'Alignment Length', 'Sequence Identity', 
                      'Query Start', 'Subject Start', 'Subject End', 'Sequence']
            writer.writerow(titles)

            for ID in output:
                for contig in output[ID]:
                    for match in output[ID][contig]:
                        row = [ID, contig] + list(match)
                        writer.writerow(row)
        else:
            # Write header
            titles = ['ID', 'Query', 'Contig', 'E-value', 'Alignment Length', 'Sequence Identity', 
                      'Query Start', 'Subject Start', 'Subject End', 'Sequence']
            writer.writerow(titles)

            for ID in output:
                for query in output[ID]:
                    for match in output[ID][query]:
                        row = [ID, query] + list(match)
                        writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse BLAST output and export to CSV.')

    parser.add_argument('--directory', type=str, required=True, help='Path to directory with BLAST .out files')
    parser.add_argument('--querys', type=int, required=True, help='Number of queries (1 or multiple)')
    parser.add_argument('--blast_type', type=str, choices=['blastn', 'tblastn'], required=True, help='BLAST type: blastn or tblastn')
    parser.add_argument('--max_Eval', type=float, default=10.0, help='Maximum E-value threshold')
    parser.add_argument('--sequence_identity', type=int, default=0, help='Minimum sequence identity percentage')
    parser.add_argument('--min_length', type=int, default=0, help='Minimum alignment length as a percentage of query length')
    parser.add_argument('--file_name', type=str, required=True, help='Output CSV file name (without extension)')

    args = parser.parse_args()

    # Ensure directory ends with "/"
    if not args.directory.endswith("/"):
        args.directory += "/"

    output = parseDir(args.directory, args.querys, args.blast_type, args.max_Eval, args.sequence_identity, args.min_length)
    output_csv(output, args.querys, args.file_name)