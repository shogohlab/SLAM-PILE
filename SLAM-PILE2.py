'''
inputs: gtf, fa, pileup output file
output: T to C conversion statistic from pileup
'''

import pandas as pd

'''gtf'''
print('\n\n\n1. gtf --- read gtf\n')

#####Denote desired gtf file.
df = pd.read_csv('desired_gtf_file.gtf', sep='\t', skiprows=5, header=None)
print(df)

'''gtf --- only select type: gene'''
print('\n\n\n2. gtf --- only select type: gene\n')
df = df[df[2]=='gene']
print(df)

'''gtf --- reduce df columns to ['ensg','chr','start','end','strand']'''
def get_ensg(string):
    return string.split(';')[0].split(' ')[1][1:-1]

print('''\n\n\n3. gtf --- reduce df columns to ['ensg','chr','start','end','strand']\n''')
ensgs = []
for info in df[8]:
    ensg = get_ensg(info)
    ensgs.append(ensg)

df['ensg'] = ensgs
df['chr'] = df[0]
df['start'] = df[3]
df['end'] = df[4]
df['strand'] = df[6]
df = df[['ensg','chr','start','end','strand']]
print(df)

'''gtf --- check gtf chromosomes'''
chrs = []
for ch in df['chr']:
    if ch not in chrs:
        chrs.append(ch)
print('\n   number of gtf chromosomes: ' + str(len(chrs)))
'''gtf --- check ensg redundancy'''
ensgs.sort()
num_ensg_redundance = 0
for i in range(1,len(ensgs)):
    if ensgs[i] == ensgs[i-1]:
        num_ensg_redundance += 1
print('\n   number of ensg redundance: ' + str(num_ensg_redundance))



'''fa, gtf --- mapping fa chrs to gtf chrs'''
print('''\n\n\n4. fa --- read fa\n''')
'''fa --- checking .fa chromosomes'''

#####Denote desired genome.fa file.
file = open('desired_genome_fa_file.fa','r')

lines = file.readlines()
chrs_fa = []
chrs_fa_sequence = []
count=0
for line in lines:
    if line[0] == '>':
        if count!=0:
            chrs_fa_sequence.append(sequence)
        ch_fa=line[1:-1]
        chrs_fa.append(ch_fa)
        count+=1
        sequence = ''
    else:
        sequence += line[:-1]
chrs_fa_sequence.append(sequence)

print('   number of fa chromosomes: ' + str(len(chrs_fa)))

print('''\n\n\n5. fa --- mapping fa chrs to fa sequences\n''')
chrs_fa_sequence_dic = {}
for i in range(len(chrs_fa)):
    chf = chrs_fa[i]
    chrs_fa_sequence_dic[chf] = chrs_fa_sequence[i].upper()

'''fa --- cross checking fa gtf chromosomes'''
print('''\n\n\n6. gtf,fa --- cross checking fa gtf chromosomes\n''')
chf_ch_dic={}
for chf in chrs_fa:
    string = chf + ': '
    if 'v' in chf:
        chf_change = chf.split('_')[1][:-2] + '.' + chf.split('_')[1][-1]
    else:
        chf_change = chf
    for ch in chrs:
        if chf_change == ch:
            string += ch
            string += ', '
            chf_ch_dic[chf] = ch

ch_chf_dic={}
for chf in chf_ch_dic.keys():
    ch = chf_ch_dic[chf]
    ch_chf_dic[ch] = chf
print('\n   getting gtf2fa & fa2gtf chromosomes dictionary')
print('   number of gtf2fa dictionary chromosomes: ' + str(len(ch_chf_dic)))
print('   number of fa2gtf dictionary chromosomes: ' + str(len(chf_ch_dic)))

'''gtf,fa --- link gtf.ensgs to fa.sequence'''
print('''\n\n\n7. gtf,fa --- link gtf.ensgs to fa.sequence\n''')

'''map ensg to gtf location in such format: chr:start-end'''
df2=df.set_index('ensg')
ensg_location_dic={}
ensg_strand_dic={}
for ensg in ensgs:
    location = df2['chr'][ensg] +':'+ str(df2['start'][ensg]) +'-'+ str(df2['end'][ensg])
    ensg_location_dic[ensg]=location
    ensg_strand_dic[ensg]=df2['strand'][ensg]

'''map gtf location to fa.sequence'''
def parse_location(location):
    chromosome, start, end = location.split(':')[0], int(location.split(':')[1].split('-')[0]), int(location.split(':')[1].split('-')[1])
    return chromosome, start, end

def get_sequence_from_ensg(ensg, ensg_location_dic, ch_chf_dic, chrs_fa_sequence_dic):
    test_location = ensg_location_dic[ensg]
    chromosome, start, end = parse_location(test_location)
    chromosome = ch_chf_dic[chromosome]#switch to fa chromosome
    sequence = chrs_fa_sequence_dic[chromosome][start-1:end]
    return sequence

def get_sequence_from_location(location, ch_chf_dic, chrs_fa_sequence_dic):
    chromosome, start, end = parse_location(location)
    chromosome = ch_chf_dic[chromosome]#switch to fa chromosome
    sequence = chrs_fa_sequence_dic[chromosome][start-1:end]
    return sequence

def get_Tcontent(sequence):
    Tcontent = 0
    for base in sequence:
        if base == 'T':
            Tcontent += 1
    return Tcontent

'''adding gene length and Tcontent to df'''
print('''\n\n\n8. adding gene length and Tcontent to df\n''')
length, Tcontent = [], []
num_keyerror = 0
for ensg in df['ensg']:
    try:
        seq = get_sequence_from_ensg(ensg, ensg_location_dic, ch_chf_dic, chrs_fa_sequence_dic)
        length.append(len(seq))
        Tcontent.append(get_Tcontent(seq))
    except KeyError:
        num_keyerror += 1
        length.append('nif')
        Tcontent.append('nif')
        
df['length'] = length
df['Tcontent'] = Tcontent
print(df.shape)
df = df[df['length']!='nif']
print(df)
print('number of key errors: '+str(num_keyerror))

'''output area'''
print('\n\n\nmajor outputs:\n')
print('1. df')
print('''   gtf as dataframe with columns: ['ensg','chr','start','end','strand','length','Tcontent']''')
print('2. ch_chf_dic')
print('   gtf chromosome ---> fa chromosome')
print('3. chrs_fa_sequence_dic')
print('   fa chromosome ---> fa sequence')
print('4. ensg_location_dic')
print('   ensg ---> fa location')
print('\n')

'''filter out non Ts'''
print('''\n\n\n9. generating coverageOnTs & conversionsOnTs\n''')

def base_checker(position, chromosome_sequence, base='T'):
    fa_base = chromosome_sequence[position-1]
    if fa_base == base:
        return True
    else:
        return False
    
def get_pileup_count(sel_df2, conversion=False, strand='+'):
    if conversion==True:
        if strand == '+':
            sel_df2 = sel_df2[sel_df2['nucleotide'] == 'C']
        elif strand == '-':
            sel_df2 = sel_df2[sel_df2['nucleotide'] == 'G']
    else:
        pass
    count = sum(sel_df2['count'])
    return count

def bisect_search(value_list, target_value, lower_index, upper_index, mode = 'find_end'):

    '''this is a recursive funciton to find the index of target_value in a sorted ascending value_list.
        
        if target_value doesn't exist in list:
            if mode = 'find_end': return the index with maximum value < target_value;
            if mode = 'find_start': return the index with minimum value > target_value;'''

    try:
        if value_list[lower_index] == target_value:
            return lower_index
        if value_list[upper_index] == target_value:
            return upper_index
        else:
            if upper_index - lower_index == 0:
                return lower_index
            if upper_index - lower_index == 1:
                if mode == 'find_end':
                    if value_list[upper_index] < target_value:
                        return upper_index
                    else:
                        return lower_index
                if mode == 'find_start':
                    if value_list[lower_index] > target_value:
                        return lower_index
                    else:
                        return upper_index
            else:
                mid_index = int((lower_index + upper_index)/2)
                if value_list[mid_index] == target_value:
                    return mid_index
                else:
                    if value_list[mid_index] > target_value:
                        upper_index = mid_index
                    else:
                        lower_index = mid_index
                    return bisect_search(value_list, target_value, lower_index, upper_index, mode = mode)
    except IndexError:
        return -1

filter_it = True

#####Denote sample filename here.
filename = 'filename'

pileup_filename = filename+'.pileup.csv'
print(pileup_filename + '\n')

if filter_it:
    pileup_df = pd.read_csv(pileup_filename)
    print(pileup_df)

    print('\nfiltering for all the Ts\n')
    current_chromosome = ''
    t_index=[]
    a_index=[]
    for i in range(pileup_df.shape[0]):
        chromosome, position = pileup_df['seqnames'][i], pileup_df['pos'][i]
        strand = pileup_df['strand'][i]
        
        if chromosome != current_chromosome:
            current_chromosome = chromosome
            chromosome_sequence = chrs_fa_sequence_dic[chromosome]

        if strand == '+':
            if base_checker(position, chromosome_sequence, base='T') == True:
                t_index.append(i)
        elif strand == '-':
            if base_checker(position, chromosome_sequence, base='A') == True:
                a_index.append(i)
        

    pileup_df_t = pileup_df.iloc[t_index]
    print(pileup_df_t)
    pileup_df_t.to_csv(pileup_filename[:-4] + '.t' + pileup_filename[-4:],index=False)

    pileup_df_a = pileup_df.iloc[a_index]
    print(pileup_df_a)
    pileup_df_a.to_csv(pileup_filename[:-4] + '.a' + pileup_filename[-4:],index=False)

else:
    pileup_df_t = pd.read_csv(pileup_filename[:-4] + '.t' + pileup_filename[-4:])
    pileup_df_a = pd.read_csv(pileup_filename[:-4] + '.a' + pileup_filename[-4:])



print('\nstart counting\n')
coverageOnTs, conversionsOnTs = [], []
current_chromosome = ''

for ensg in df['ensg']:
    print(ensg)
    try:
        location = ensg_location_dic[ensg]
        strand = ensg_strand_dic[ensg]
        chromosome, start, end = parse_location(location)
        chromosome = ch_chf_dic[chromosome]
        if chromosome != current_chromosome:
            current_chromosome = chromosome
            if strand == '+':
                sel_df1 = pileup_df_t[pileup_df_t['seqnames'] == chromosome]
            elif strand == '-':
                sel_df1 = pileup_df_a[pileup_df_a['seqnames'] == chromosome]
            value_list = list(sel_df1['pos'])
            lower_index = 0
            sel_df1_row_num = sel_df1.shape[0]
            upper_index = sel_df1_row_num - 1

        target_value_start = start
        target_index_start = bisect_search(value_list, target_value_start, lower_index, upper_index, mode = 'find_start')
        target_value_end = end
        try:
            target_index_end = bisect_search(value_list, target_value_end, target_index_start, upper_index, mode = 'find_end')
        except RecursionError:
            print(sel_df1.iloc[target_index_start-5])
            print(sel_df1.iloc[target_index_start-3])
            print(sel_df1.iloc[target_index_start-2])
            print(sel_df1.iloc[target_index_start-1])
            print(sel_df1.iloc[target_index_start])
            print(sel_df1.iloc[target_index_start+1])
            print(sel_df1.iloc[target_index_start+2])
            print(sel_df1.iloc[target_index_start+3])
            print(sel_df1.iloc[target_index_start+4])
            break

        sel_df2 = sel_df1.iloc[target_index_start : target_index_end+1]
        coverage = get_pileup_count(sel_df2, conversion=False, strand=strand)
        conversion = get_pileup_count(sel_df2, conversion=True, strand=strand)

        coverageOnTs.append(coverage)
        conversionsOnTs.append(conversion)

    except KeyError:
        coverage = 'nif'
        conversion = 'nif'
        coverageOnTs.append(coverage)
        conversionsOnTs.append(conversion)

    print(coverage)
    print(conversion)

df['coverageOnTs'] = coverageOnTs
df['conversionsOnTs'] = conversionsOnTs
conversionRate = []
for i in range(len(coverageOnTs)):
    coverage = coverageOnTs[i]
    if coverage == 0:
        conversionRate.append('nil')
    else:
        conversion = conversionsOnTs[i]
        cr = conversion/coverage
        conversionRate.append(cr)
df['conversionRate'] = conversionRate

print(df)

df.to_csv(pileup_filename[:-4] + '.conversionStatistic' + pileup_filename[-4:],index=False)




