# run this script before running rxn.py with the option
# exclude_species_without_temp_data = True

types = ['sec_species','gas_species','minerals']

fh = open('hanford.dat','r')
h_lines = fh.readlines()
fh.close()

num_temps = 8
deg_25 = 1

blocks = []
block = []
for i,line in enumerate(h_lines):
    if line.split()[0] == "'null'":
        blocks.append(block)
        block = []
    else:
        block.append(line)

for i,t in enumerate(types):
    print t
    lines = blocks[i+1]
    fh = open(t + '_no_temp_data.txt','w')

    for line in lines:
        f = line.strip().split()
        name = f[0]
        if 'sec' in t:
            idx = 1
        else:
            idx = 2
        num = int(f[idx])
        temp_data = [float(x.replace('d','e').replace('D','E')) for x in f[2*num+idx+1:2*num+num_temps+1+idx]]

        non25_temp_data = []
        non25_temp_data.extend(temp_data[0:deg_25])
        non25_temp_data.extend(temp_data[deg_25+1:])
        if sum([abs(x-500.0) for x in non25_temp_data]) < 1.0:
            fh.write(name +'\n')
            
    fh.close()
