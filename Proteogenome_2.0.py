#!usr/bin/python3

import numpy as np
import pandas as pd
import re

class HCMV:
          
    def __init__(self, FASTA_filename='', annot_filename='', upload_files=True):
        """
        Version: 1.1
        Author: Giammarco Ferrari
        Thank you for using Proteogenome!!!
        """
        if FASTA_filename: self.FASTA_filename=FASTA_filename
        if annot_filename: self.annot_filename=annot_filename
        self.prot_index = {} # Prepare the data structure for the protein index

        if upload_files:
            self.FASTA_lst=self.file_to_lst(self.FASTA_filename)
            self.annot_lst=self.file_to_lst(self.annot_filename)
    	

    def file_to_lst(self, file_name):
        """
        Version : 1.0 

        Name History : file2list (from Proteogenome 1.0)

        This function upload in a list the file rows.
        For each row it is removed the last character that usually is a ('/n'). 

        INPUT  :  file_name   String   The name of the file to upload.
        
        OUTPUT :  file_in_lst List     The output list with the content of the file.
        """
        file_in_lst = []
        file_hand = open(file_name,'r')

        for row in file_hand:
            if len(row) != 0:
                file_in_lst.append(row[:-1])

        file_in_lst=file_in_lst[:-1] # I do not consider the last row of the file because the previous
                                     # operation row[:-1] generates an empty string.
        return file_in_lst

    
    def print_lst(self, input_list, limit):
        """
        Version : 1.0

        Print the list content limitely to the element indicated by the parameter 'limit'. 
        INPUT  : input_list     List      Items to print
                 limit          Int       Number of items to print 
        
        OUTPUT :
        """
        for ind, i in enumerate(input_list):
            if ind < limit:
                print(i)

    def annot_to_df(self, annotations, annot_format ='gff3'):
        """
        Version: 1.0
        Name History: annot_to_df
        Extract relevant data from the annotation file and put them in a DataFrame
        INPUT : annotations     List    List list filled with the annotation rows
        OUTPUT:
        """
        #annotations_fh = open(GFF3,'r')

        protein_index = {}
        self.annotations_df = pd.DataFrame()

        new_row = {}      # In order to append a new row to the dataframe I need a dictionary 

        # **************** Parsing file for GFF3 format
        if annot_format =='gff3':
            unique_pat = re.compile(r'(.*?)\t')
            ID_pat = re.compile(r'.*?\tID=(.*?);')            # specific patterns for the annotation of the GFF3 format
            gene_pat = re.compile(r'.*?;gene=(.*?);')
            product_pat = re.compile(r'.*?;product=(.*?)$')
            col_index=[2, 3, 4, 6]

        elif annot_format == 'gtf':
            unique_pat = re.compile(r'(.*?)\t')
            pass


        #for row in annotations_fh:
        for row in annotations:

            if (row[0]!='#') and ('country=United Kingdom: Cardiff;culture-collection=ATCC' not in row): 
                match_lst=unique_pat.findall(row)
                
                # print(row)
                # a=input()
                try:
                    new_row['row_type'] = match_lst[col_index[0]]      # 2 - GFF3 column 3
                    new_row['coordinate1'] = match_lst[col_index[1]]   # 3 - GFF3 column 4
                    new_row['coordinate2'] = match_lst[col_index[2]]   # 4 - GFF3 column 5
                    new_row['strand'] = match_lst[col_index[3]]        # 6 - GFF3 column 7
                except:
                    print('This line returns an error on parsing procedure!!!')
                    print(f'-------------------------------------------------\n {row}')
                    print('Please press a button to continue')
                    #a=input()
                
                try:
                    new_row['ID'] = ID_pat.match(row).group(1)
                except:
                    new_row['ID'] = 'ID PASS'
                try:
                    new_row['gene'] = gene_pat.match(row).group(1)
                except:
                    new_row['gene'] = 'gene PASS'
                try:
                    new_row['product'] = product_pat.match(row).group(1)
                    new_row['product'] = new_row['product'].split(';')[0]
                except:
                    new_row['product'] = 'product PASS'  

                #print(new_row) 

                self.annotations_df = self.annotations_df.append(new_row, ignore_index=True)
        self.annotations_df = self.annotations_df[['row_type', 'coordinate1', 'coordinate2', 
                                                   'strand', 'ID', 'gene', 'product']]
        #annotations_fh.close()

        #return self.GFF3_to_df
       


    def prot_index(self, annot_df):
        """
        Version:
        Name History:
        Generate a dictionary of proteins ID ('gene' field) as key and the exon structure in a list as value.
            Example:
             {'RL1'   : []
              'RL5A'  : []
              'RNA2.7': []
              .....
              }

        INPUT :
        OUTPUT:
        """
        #self.prot_index = {}

        #print(annot_df['gene'])
        prot_gene_lst=set(annot_df['gene'].tolist())
        # print(len(prot_gene_lst))

        
        for prot_gene in prot_gene_lst:
            #if prot_gene=='UL150A':print(f'--- --- {prot_gene} BLOCK --- --- ')
            #print(prot_gene)
            if prot_gene != 'gene PASS': # Only if exist a gene (otherwise means that it is a non-coding region)
                # Extract the exons that belong to the gene structure 
                #print('extract prot_block')
                prot_block = annot_df[(annot_df['gene']==prot_gene) & (annot_df['row_type']=='CDS')][['coordinate1', 'coordinate2', 'strand']]
                #print('prot_block --- ',prot_block)
                for exon in prot_block.iterrows(): # Coordinates and strand of the exons in the each protein come from a groupby object 
                    #if prot_gene=='UL150A':print('exon --- ',exon)
                    #print('iterate prot_block')
                    exon_feat=[]                   # What I want is ti have a list of lists where avery sublist is a exon representation 
                    for feat in exon[1]:   # exon is a pandas serie and the index 1 get the table of the exons           
                        exon_feat.append(feat)
                        #print('append exon_feat')
                    #if prot_gene=='UL150A':print('exon_feat --- ',exon_feat)
                    if (prot_gene in self.prot_index):
                        if (exon_feat!=[]) & (exon_feat not in self.prot_index[prot_gene]):
                            self.prot_index[prot_gene].append(exon_feat) # Append the features of the actual exon on the list 
                                                                         # of exons that belongs to the actual protein.
                    else:
                        self.prot_index[prot_gene]=[exon_feat]

                #if prot_gene=='UL150A':print(exon_feat)
                #a=input()

            #return self.prot_index


# ****************************************************************************************** #
# ******************************* FILE MANIPULATION FOR PoGo ******************************* #

    
    def rectify_rows(self, list_rows, target_sub_str=[], target_patterns=[]):
        """
        This function recive a list of rows. Therefore cleans every single row replacing the 
        target substrings with the proper substitution.
        Target substrings can be substitute in two ways:
                - precise patterns: are substrings known before the substituting operation 
                                    and are applied through the command ".replace()".
                - generic patterns: substrings that can change in the but with the same 
                                    repeated structure. This type of substring must be set as
                                    regerx Python patterns enclosed in ''.
                                      
        INPUT : list_rows           List    List of rows to clean
                target_sub_str      List    List of tuples of two elements. Replacement through 
                                            .replace() command. 
                                            Tuple position reference:
                                            [0] target substring that must be substitute.
                                            [1] replacement string used in the substitution
                target_patterns     List    List of tuples of two elements. Replacement through 
                                            .sub() command in a regex search.
                                            Tuple position reference:
                                            [0] target pattern that must be substitute.
                                            [1] replacement string used in the substitution
        OUTPUT: list_rows           List    The original list of rows cleaned from the 
                                            substrings provided as input.
        """
        
        for ind, row in enumerate(list_rows):
            #print(row)
            if target_sub_str:
                for clear_this in target_sub_str:
                    row = row.replace(clear_this[0],clear_this[1])

            
            if target_patterns:
                for clear_this in target_patterns:
                            #        |target pattern|    |replacement sub|  
                    row = re.sub(r'' + clear_this[0] + '', clear_this[1], row) 
            
            #print(row)
            #a=input()    
            list_rows[ind] = row
            #print('\nlist_rows\n',list_rows)

        return list_rows
    


    def check_format(self, list_rows):
        """
        Version: 1.0

        Name History: check_format

        This function controls if the first column in a datafile contains the expected data.
        The default is the GFF3 control where '##' chars represent comment lines and 
        the first column in a .split() row contains the strain identifier.

        INPUT :
        OUTPUT:
        """
        format_control=True    # I suppose that the file has the expected format
        strain_ID='NA'

        for row in list_rows:
            #print(row[:2])             # Find the strain ID as future refernce 
            if row[:2] != '##':           # The '##' represent comment lines in the GFF3 format
                strain_ID=row.split()[0]  # In a normal row the strain ID is in the first column
                break
        print(f'Row header - {strain_ID}')
        for row in list_rows:             # Restart from the beginning of the list 
            
            if row[:2] != '##':
                if row:
                    current_strain_ID=row.split()[0]
                    if current_strain_ID != strain_ID:
                        print('INVALID ROW')
                        print(row)
                        format_control=False
                        break
        if format_control:
            print('FILE FORMAT - OK')
        else:
            print('XXXXXX  WARNING- WRONG FILE FORMAT  XXXXXX') 

        print('++ FILE FORMAT CONTROL FINISHED ++')


    def locus_tag_substitution(self, FASTA_lst):
        """
        Version: 1.0

        Name History: locus_tag_substitution

        This function substitutes the tag 'locus_tag='' with 'gene:gene-' and 'transcript:rna-'
        into the FASTA rows. 
        This function is specific to the single operation of fitting the FASTA header into 
        a PoGo readable format.

        INPUT :     FASTA_lst       List    List of FASTA rows including also FASTA sequences. 
        OUTPUT:     FASTA_lst4PoGo  List    List of FASTA rows with the rectified headers and
                                            the FASTA sequences.
        """
        FASTA_lst4PoGo=[]

        locus_tag_patt=re.compile(r'.*?locus_tag=(.*?)\s')
        for row in FASTA_lst:
            if row[0] =='>':
                target_tag=locus_tag_patt.search(row).group(1)
                #print(row)
                row=row.replace('locus_tag='+target_tag, 
                                'gene:'+target_tag+' transcript:'+target_tag)
                                #gene:gene-   The first suggestion was to add gene:gene-
                                #             But later       
            FASTA_lst4PoGo.append(row)
                #print(row)
                #a=input()
        return FASTA_lst4PoGo    



    def protein_set_bed(self, prot_list=[], bed_fn='protein_track.bed'):
        """
        Receive a list of protein codes and create the .bed track with the genomic locations of the protins.
        Example of multi exon protein in input:
          UL37 - [['52573', '53060', '-'], ['51302', '51344', '-'], ['50262', '51197', '-']]
        INPUT:
        OUTPUT:
        """
        
        if prot_list==[]:
            prot_list=self.prot_index.keys()
            #prot_list=set(prot_list)

        print(f'Start processing {len(prot_list)} proteins')
        prot_track_fh = open(bed_fn,'w')

        chr_name = 'NC_006273.2'

        BED_rows = []
        prot_row = ''

        for protein in prot_list:
            
            prot_row=''

            chromStart = self.prot_index[protein][-1][0]    # Set the START in chromosome at the lowest coordinate in the current gene structure
            chromEnd = self.prot_index[protein][0][1]      # Set the END in chromosome at the higest coordinate in the current gene structure
            
            tickStart = self.prot_index[protein][-1][0]
            tickEnd = self.prot_index[protein][0][1]

            score = '1000'
            
            itemRGB = "255.0.0"
            blockCount = str(len(self.prot_index[protein]))
            
            blk_size_int=[]
            blockSizes_str = ''
            blockStarts_str = ''
            for exon in reversed(self.prot_index[protein]): # Iterate over the list of lists of exons from the neares exon toward the farest
                strand = exon[2]
                
                # if not chromStart: 
                #     chromStart = int(exon[0])   # If chromeStart is == None means that this is the first exon considered
                # elif chromStart > int(exon[0]): # Otherwise I have to chose the lower coordinate in the exon structure
                #     chromStart = int(exon[0])

                # if not chromEnd: 
                #     chromEnd = int(exon[1])     # If chromeEnd is == None means that this is the first exon considered
                # elif chromEnd > int(exon[1]):   # Otherwise I have to chose the higer coordinate in the exon structure
                #     chromEnd = int(exon[1])

                blk_s = int(exon[0])
                blk_e = int(exon[1])
                blk_size = str(abs(blk_s-blk_e))
                blockSizes_str += blk_size + ','

                blockStarts_str += str(abs(int(chromStart)-blk_s)) + ','
 
            blockSizes_str = blockSizes_str[:-1] # Remove the exceded ','
            blockStarts_str = blockStarts_str[:-1]

            chromStart = str(chromStart)
            chromEnd = str(chromEnd)

            prot_row += (chr_name+'\t'+chromStart+'\t'+chromEnd+'\t'+protein+'\t'+score+'\t'+
                         strand+'\t'+tickStart+'\t'+tickEnd+'\t'+itemRGB+'\t'+blockCount+'\t'+
                         blockSizes_str+'\t'+blockStarts_str+'\n')

            BED_rows.append(prot_row)

        for row in BED_rows:
            prot_track_fh.write(row)
        prot_track_fh.close()

            #print(f'{protein} - {self.prot_index[protein]}')


# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/ #
#                                              M - A - I - N                                           #
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/ #

if __name__ == '__main__':
    
    GFF3='C:/Users/Giammarco Ferrari/Documents/PoGo4Virus/Source HCMV files/HCMV_ComplRecord.gff3'
    FASTA='C:/Users/Giammarco Ferrari/Documents/PoGo4Virus/Source HCMV files/HCMV_CodingSeq_AminoAcidsFASTA.txt'
    
    print("""\n\n
                   *****************
                   * LOADING FILES *
                   *****************       
          """)


    print(GFF3)
    HCMV_instance=HCMV(FASTA, GFF3)
    print(f'String representation of the object ==> {HCMV_instance}')
    print(f'Filename for the FASTA file         ==> {HCMV_instance.FASTA_filename}')
    print(f'Filename for the GFF3  file         ==> {HCMV_instance.annot_filename}\n')

    print("""\n\n
                   *********************
                   * FASTA PREPARATION *
                   *********************       
          """)

    print('***** FASTA sequence in LIST *****\n')
    FASTA_lst=HCMV_instance.FASTA_lst
    HCMV_instance.print_lst(FASTA_lst, 30)
    print('\n' * 2)

    print('+++++++++ FASTA rectification (PoGo COMPLIANT) +++++++++\n')
    FASTA_lst=HCMV_instance.rectify_rows(FASTA_lst, target_sub_str=[('lcl|NC_006273.2_prot_',''),('[',''),
                                                                    (']',''),('..',':'),
                                                                    ('gene-',''), ('rna-','')]) 
    FASTA_lst=HCMV_instance.locus_tag_substitution(FASTA_lst)
    
    HCMV_instance.print_lst(FASTA_lst, 30)

    print("""\n\n
                    **************************
                    * ANNOTATION PREPARATION *
                    **************************       
          """)
    
    print('***** ANNOTATIONS GFF3 in list *****\n')
    GFF3_lst=HCMV_instance.annot_lst
    HCMV_instance.print_lst(GFF3_lst, 30)
    print('\n' * 2)

    HCMV_instance.check_format(GFF3_lst)

    print("""\n\n
                   ***********************
                   * USE AGAT TO CONVERT *
                   *    GFF3 INTO GTF    *
                   ***********************       
          \n\n""")

    GTF='C:/Users/Giammarco Ferrari/Documents/PoGo4Virus/HHV5_ge_ex_CD_mR_Gid_corrgene.gtf'

    print('***** ANNOTATIONS GTF in list *****\n')
    GTF_lst=HCMV_instance.file_to_lst(GTF)

    HCMV_instance.check_format(GTF_lst)

    GTF_lst=HCMV_instance.rectify_rows(GTF_lst, target_sub_str=[('gene-',''), ('rna-','')])

    HCMV_instance.print_lst(GTF_lst, 30)

    print(f'Uploading the GFF3 file in a df')
    HCMV_instance.annot_to_df(GFF3_lst)
    print(HCMV_instance.annotations_df)
    print(f'**done** \n')

    print(f'Creating Protein index')
    print(' ************ DISABLED IN THE CODE ************ ')
    # #HCMV_instance.prot_index(HCMV_instance.GFF3_to_df)
    # print(HCMV_instance.prot_index)
    # print(f'**done** \n')

    print('Create .bed track for the whole protein set of the HCMV')

    # protein_list_test=['RL1', 'UL29', 'UL150A']
    # HCMV_instance.protein_set_bed()
    # print('hello world ....... at the end!!!!!')

