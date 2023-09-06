import sys
def fix_cif_file(file, output):
    # Read file by lines
    file = open(file, 'r')
    output = open(output,'w')
    lines = file.readlines()
    # Set up dummy variable
    previous = "0"
    # Set up count to get the column order
    count= -1
    # Iterate through the lines
    for line in lines:
        # Add counts each line
        count += 1
        # Reset counts at the start of a loop_, to get the column order
        if "loop_" in line:
            count = -1
        # Get the order for each attribute
        if "esd" not in line[-4:].lower():
            if "_atom_site.auth_seq_id" in line.lower():
                seqid = count
        if "ATOM" in line[:10] or "HETATM" in line[:10]:
            line_split = line.split()
            seqid_value = line_split[seqid]
            if seqid_value.isdigit():
                previous = int(seqid_value)
                output.write(line)
            else:
                previous = str(int(previous)+1)
                new_line = line.replace(seqid_value, previous)
                output.write(new_line)
        else:
            output.write(line)
            continue

fix_cif_file(sys.argv[1], sys.argv[2])