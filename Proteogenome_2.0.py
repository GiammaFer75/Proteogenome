#!usr/bin/python3

# Proteogenome
# Version: 2.3
#
# Author: Giammarco Ferrari
#
# Last Update: 27/04/2022
# it's time to package https://www.tutorialsteacher.com/python/python-package

import numpy as np
import pandas as pd
from random import randrange
from matplotlib import colors
import re



class Organism:
          
    def __init__(self, FASTA_filename='', annot_filename='', 
                 input_table_filename='', upload_files=True):
        # """
        # Version: 1.1
        # Author: Giammarco Ferrari
        # Thank you for using Proteogenome!!!
        # """
        # 
        if FASTA_filename!='': 
            self.FASTA_filename=FASTA_filename
            self.FASTA_lst=self.file_to_lst(self.FASTA_filename)
        if annot_filename!='': 
            self.annot_filename=annot_filename
            self.annot_lst=self.file_to_lst(self.annot_filename)
        if input_table_filename!='':
            self.input_table_filename=input_table_filename
            self.input_table=self.load_input_table(self.input_table_filename)

        self.prot_CDS_index = {}      # Initialise dictionary for protein ---> CDS index
        self.prot_pep_index = {}      # Initialise dictionary for protein ---> peptide index
        self.pep_prot_index = {}      # Initialise dictionary for peptide ---> protein index
        self.prot_PSMint_index = {}   # Initialise dictionary for protein ---> PSM - intensity - RGB intensity index     
    	

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

    
    def print_lst(self, input_list, limit=0, en_sep=True, sep_type='-'):
        """
        Version : 2.0

        Name History : print_lst (from Proteogenome 1.0)

        Print the list content limitely to the element indicated by the parameter 'limit'. 

        INPUT:  input_list  List    Items to print.
                limit       Int     Number of items to print. 
                en_sept     Bool    Enable the separation of each element of the list with a string.
                sep_type    Char    The character that will make up the separator string.   

        
        OUTPUT :
        """

        for ind, i in enumerate(input_list):
            if ind < limit:
                print(i)
                if en_sep: print(sep_type * 100)

    def load_input_table(self, filename, sep='\t'):
        """
        Version: 1.0

        Name History: load_input_table

        INPUT :
        OUTPUT:
        """
        
        input_tab=np.array([['','','','','']])

        fh=open(filename, 'r')
        for row in fh:
            row=row.split(sep)
            row=np.array(row)
            row[-1]=row[-1].replace('\n','')
            input_tab=np.concatenate([input_tab,[row]], axis=0)
        fh.close()
        input_tab=input_tab[1:,:]   # Remove initialisation row
        return input_tab


    def load_generic_table(self, filename, sep='\t'):
        """
        Version: 1.0

        Name History: load_generic_table

        INPUT :
        OUTPUT:
        """
        fh=open(filename, 'r')
        first_line=fh.readline().rstrip()             # Divide in columns the first row 
        number_of_columns=len(first_line.split(sep))  # Fetch the number of columns from the first row 
        
        fh.close()

        input_tab=np.empty((1,number_of_columns))

        fh=open(filename, 'r')
        for row in fh:
            row=row.split(sep)
            row=np.array(row)
            row[-1]=row[-1].replace('\n','')
            input_tab=np.concatenate([input_tab,[row]], axis=0)
        fh.close()
        input_tab=input_tab[1:,:]   # Remove initialisation row
        return input_tab

    def annot_to_df(self, annotations, annot_format ='gff3'):
        """
        Version: 1.0

        Name History: annot_to_df
        
        This function extracts relevant data from the annotation file and put them in a DataFrame
        
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
            ID_pat = re.compile(r'.*?\tID=(.*?);')            # specific patterns for the 
            gene_pat = re.compile(r'.*?;gene=(.*?);')         # annotation in GFF3 format
            product_pat = re.compile(r'.*?;product=(.*?)$')
            col_index=[2, 3, 4, 6]

        elif annot_format == 'gtf':
            unique_pat = re.compile(r'(.*?)\t')
            ID_pat = re.compile(r'.*?\sgene_id\s\"(.*?)\";')   # specific patterns for the 
            gene_pat = re.compile(r'.*?;\sgene\s\"(.*?)\";')   # annotation in GTF format
            product_pat = re.compile(r'.*?;\sproduct\s\"(.*?)\"$')
            col_index=[2, 3, 4, 6]

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
                    new_row['ID'] = 'ID nf'
                try:
                    new_row['gene'] = gene_pat.match(row).group(1)
                except:
                    new_row['gene'] = 'gene nf'
                try:
                    new_row['product'] = product_pat.match(row).group(1)
                    new_row['product'] = new_row['product'].split(';')[0]
                except:
                    new_row['product'] = 'product nf'  

                #print(new_row) 

                self.annotations_df = self.annotations_df.append(new_row, ignore_index=True)
        self.annotations_df = self.annotations_df[['row_type', 'coordinate1', 'coordinate2', 
                                                   'strand', 'ID', 'gene', 'product']]
        #annotations_fh.close()

        #return self.GFF3_to_df
       
                

    def CDS_annot_matrix(self, annotations):
        """
        Version: 1.0

        Name History: annot_matrix, CDS_annot_matrix

        This function receives a list of annotations rows and extracts only the block of
        rows regarding the protein CDS 
        INPUT :
        OUTPUT:
                self.annotation_matrix
                self.CDS_matrix
        """

        first=True
        for row in annotations:
            if row[0]!='#':
                                                     
                if first:
                    row_array=np.array(row.split('\t'),dtype=object)          # Split the current annotation rows in columns
                    self.annotation_matrix = np.array([row_array]) # increase the annotation matrix dimension
                    first=False
                else:
                    row_array=np.array(row.split('\t'),dtype=object)
                    
                    self.annotation_matrix = np.concatenate((self.annotation_matrix, [row_array]))
                    # print('*'*50)
                    # print(self.annotation_matrix)
                    # a=input()

        self.CDS_matrix=self.annotation_matrix[self.annotation_matrix[:,2]=='CDS'] # Filter only for coding regions 


    def protein_CDS_index(self, annotations, annot_format='gff3'):
        """
        Version: 3.0

        Name History: prot_index, protein_CDS_index
        
        Generate a dictionary with protein ID ('gene' field in GFF3) as key and 
        its CDS list as the value.

            Example:
             {'RL1'   : [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
              'RL5A'  : [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
              'RNA2.7': [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
              .....
              }
        
        Two alternative patterns could be used to parse the protein identifier 
        into the annotation file:
            - GFF3 ----> .*?;gene=(.*?);
            - GTF  ----> .*?;\sgene\s\"(.*?)\";
          

        INPUT : annotations     [List]  The annotation file stored in a list
                annot_format    [Str]   Tag that indicates which kind of format the annotation file is.

        OUTPUT: .prot_index     [Dict]  The protein index with the CDS coordinates.
                                        (Key)   [Str]               Protein ID
                                        (Value) [List][List][Str]   CDS coordinates = 
                                                                    chr,start,end,strand  
        """
        #self.prot_CDS_index = {}
        
        # **************** Parsing annotations in GFF3 format
        if annot_format =='gff3':                             # specific patterns for the 
            gene_pat = re.compile(r'.*?;gene=(.*?);')         # annotation in GFF3 format


        # **************** Parsing annotations in GTF format
        elif annot_format == 'gtf':                           # specific patterns for the 
            gene_pat = re.compile(r'.*?;\sgene\s\"(.*?)\";')   # annotation in GTF format

        for row in annotations:

            strand=row[6]
            
            if strand=='+':
                coord_1=row[3]
                coord_2=row[4]
            else:             # If the coding region is in reverse starnd
                coord_1=row[4]
                coord_2=row[3]

            protein_ID=gene_pat.match(row[-1]).group(1) # Find the current protein identifier.
            # print(row)
            # print(protein_ID)
            # print('-----------------------')   

            CDS_feat=[coord_1, coord_2, strand]

            if (protein_ID in self.prot_CDS_index):
                # Append the features of the current CDS on the list of CDS that belongs to the current protein.
                self.prot_CDS_index[protein_ID].append(CDS_feat) 
                                                                  
            else:
                self.prot_CDS_index[protein_ID]=[CDS_feat]

    
    def groupby_protein(self, pep_input_table, protein_ID):
        """
        Version: 1.0

        Name History: groupby_protein

        INPUT :
        OUTPUT:
        """
        protein_peptide_block=pep_input_table[pep_input_table[:,0]==protein_ID]
        return protein_peptide_block

    def protein_peptide_index(self, pep_input_table): 
        """
        Version: 1.0

        History Name: protein_peptide_index

        This function parse the input table and group all the peptides by protein ID

        INPUT : pep_input_table         np.array    
        
        OUTPUT: self.prot_pep_index     [Dict]  The protein index with the CDS coordinates.
                                            (Key)   [Str]               Protein ID
                                            (Value) [List][List][Str]   peptides = 
                                                                        pepsequence, intensisty
        """
        
        prot_ID_array = np.unique(pep_input_table[:, 0]) # Extract the protein IDs
        print(prot_ID_array)
         
        for protein in prot_ID_array:
            peptides_block =  self.groupby_protein(pep_input_table, protein) # Group the peptides that belong to the current protein.
            self.prot_pep_index[protein] = peptides_block

    
    def groupby_peptide(self, pep_input_table, peptide_sequence):
        """
        Version: 1.0

        Name History: groupby_peptide

        INPUT : pep_input_table
                peptide_sequence
        OUTPUT: protein_peptide_block
        """
        peptide_protein_block=pep_input_table[pep_input_table[:,1]==peptide_sequence]
        peptide_protein_block=peptide_protein_block[:,0]
        return peptide_protein_block

    
    def peptide_protein_index(self, pep_input_table):
        """
        Version: 1.0

        History Name: peptide_protein_index

        This function parse the input table and group all the protein ID by peptide sequence.

        INPUT : pep_input_table         np.array    
        
        OUTPUT: self.pep_prot_index     [Dict]  The peptide sequence.
                                                (Key)   [Str]       peptide sequence
                                                (Value) [List][Str] Protein IDs where the peptide come from. 
 
        """
        peptide_seq_array = np.unique(pep_input_table[:, 1]) # Extract the protein IDs
        print(peptide_seq_array)
         
        for peptide in peptide_seq_array:
            protein_block =  self.groupby_peptide(pep_input_table, peptide) # Group the peptides that belong to the current protein.
            self.pep_prot_index[peptide] = protein_block



    def protein_PSM_int_index(self):
        """
        Version: 1.0

        Name History: protein_PSM_int_index

        This function creates the protein index for the PSM and intensities. 
        Moreover generates the RGB code for each protein intensity. 
        The RGB codes will be used for the creation of the protein map.  
        
        INPUT :
                self.prot_pep_index
                self.prot_CDS_index
        OUTPUT:
        """
        import math 

        def generate_color_gradient(color_lst, reverse_gradient=False):
            """
            Version: 1.0

            Name History: generate_color_gradient
            
            ** MAIN FUNCTION ** ---> protein_PSM_int_index

            INPUT  :  color_lst     List    List of strings with the colour names that will be compose the colour gradient
                                            Example: ['black','gray','blue','green','yellow','red']
            OUTPUT :
            """

            def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
                """
                Version: 1.0

                Name History: colorFader

                ** MAIN FUNCTION ** ---> generate_color_gradient

                This function generates an array that contains an RGB gradient from color1 to the color2.
                The gradient is expressed in RGB codes.
                INPUT :
                OUTPUT:
                """
                c1=np.array(colors.to_rgb(c1))
                c2=np.array(colors.to_rgb(c2))
                return colors.to_rgb((1-mix)*c1 + mix*c2)
                

            def rgb_block(color1,color2, n_data_point):
                col = []
                for data_point in range(n_data_point): 
                    col.append(colorFader(color1,color2,data_point/n))
                return col

            n=80
            color_blocks_lst = []
            for ind, color in enumerate(color_lst):
                color_blocks_lst.append([color, color_lst[ind+1]])
                if ind == (len(color_lst)-2): break
            gradient =[]
            for block in color_blocks_lst:
                gradient += rgb_block(block[0],block[1],n)
            if reverse_gradient: gradient.reverse()
            return gradient

        def exprlev_resc_RGB(values,RGB_scale):
            """
            This function find the RGB code that 
            """
            old_max = min(values)
            old_min = max(values)
            new_max = len(RGB_scale)
            new_min = 0

            rescaled_values = []

            for old_value in values:
                try:
                    NewValue = int((((old_value - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min)
                except:
                    print('{} - {} - {}'.format(old_max, old_min, new_min))
                    a=input()
                rescaled_values.append(NewValue)

            return rescaled_values

        def vectorise_RGB_tuples(RGB_tuples, prot_exp_resc):
            
            prot_expressions_RGB=[]

            for RGB_pos in prot_exp_resc:
                RGB_toup = [round(num,3) for num in RGB_tuples[RGB_pos-1]]
                RGB_code = str(RGB_toup[0]) + ',' + str(RGB_toup[1]) + ',' + str(RGB_toup[1])
                prot_expressions_RGB.append(RGB_code)

            return prot_expressions_RGB


        # ------- MAIN ------- protein_PSM_int_index
        
        intensities=[]

        max_intensity=0
        min_intensity=0
        
        for protein, pep_array in self.prot_pep_index.items():
            PSM_sum=0
            inten_sum=0
            for pep_row in pep_array:
                PSM_sum+=int(pep_row[3])
                inten_sum+=int(pep_row[4])
            self.prot_PSMint_index[protein]=[PSM_sum, inten_sum]
            intensities.append(inten_sum)

            if max_intensity<inten_sum: max_intensity=inten_sum    
            if min_intensity>inten_sum: min_intensity=inten_sum
        
        print(f'{len(self.prot_pep_index)} - {len(self.prot_CDS_index)}')
        RGB_tup = generate_color_gradient(color_lst=['gray','blue','green','yellow','red'],
                                             reverse_gradient=False)
        prot_expressions_rescaled = exprlev_resc_RGB(intensities,RGB_tup)
        #print(len(prot_expressions_rescaled))
        prot_expressions_RGB=vectorise_RGB_tuples(RGB_tup, prot_expressions_rescaled)

        print(f'{len(prot_expressions_RGB)} - {len(self.prot_CDS_index)}')
        #print(prot_expressions_RGB)
        #print(self.prot_CDS_index)

        ind=0
        for prot, PSMint in self.prot_pep_index.items(): # Update the dictionary of self.prot_CDS_index with the RGB intensities
            self.prot_PSMint_index[prot].append(prot_expressions_RGB[ind])
            ind+=1



    def initialise_indexes(self, annot_format='gff3'):
        """
        """
        print("""
               **********************
               INDEXES INITIALISATION
               **********************
               """)
        self.CDS_annot_matrix(self.annot_lst)
        self.protein_CDS_index(self.CDS_matrix)
        self.protein_peptide_index(self.input_table)
        self.peptide_protein_index(self.input_table)
        self.protein_PSM_int_index()


# ****************************************************************************************** #
# ******************************* FILE MANIPULATION FOR PoGo ******************************* #


    def rectify_rows(self, list_rows, target_sub_str=[], target_patterns=[], 
                     open_patterns=[], rows_not_modif=False):
        """
        Version: 1.0

        Name History: rectify_rows

        This function receives a list of rows. Therefore cleans every single row replacing the 
        target substrings with the proper substitution.
        Target substrings can be substitute in two ways:
                
                - target_sub_str:   are substrings known before the substituting operation 
                                    and are applied through the command ".replace()".
                                    These patterns must be passed through the list
                                    "target_sub_string"
                
                - target_patterns:  substring that must be replaced in the file rows occur in   
                                    the same pattern. This pattern must be passed as a regerx 
                                    Python patterns enclosed in ' '. This option search for the
                                    pattern in the row. If it is found, the pattern is replaced 
                                    by the substring passed in the tuple with the pattern itself.  
                                    
                - open patterns:    In this case also the substring that will replace the target
                                    pattern is not specified explicitely as second argument in the 
                                    tuple. Instead, it must be extracted by a generic pattern as well. 

                                      
        INPUT : list_rows           [List]  List of rows to clean
                target_sub_str      [List]  List of [tuples] of two elements. Replacement through 
                                            .replace() command. 
                                            Tuple position reference:
                                            [0] target substring that must be substitute.
                                            [1] replacement string used in the substitution.
                target_patterns     [List]  List of [tuples] of two elements. Replacement through 
                                            .sub() command in a regex search.
                                            Tuple position reference:
                                            [0] target pattern that must be substitute.
                                            [1] replacement string used in the substitution.
                open_patterns       [List]  List of [tuples] of two elements. Replacement through 
                                            .sub() command in a regex search.
                                            Tuple position reference:
                                            [0] target pattern that must be substitute.
                                            [1] replacement pattern used to find the substring
                                                that will substitute the value of the pattern [0].

        OUTPUT: list_rows           [List]  The original list of rows cleaned from the 
                                            substrings provided as input.
        """ 
        rows_not_modif_lst=[]
        rows_not_modif_flag=False

        for ind, row in enumerate(list_rows):
            #print('*'*30,' NEW ROW ','*'*30)
            #print('original row - ', row)
            if target_sub_str:
                for substitution in target_sub_str:
                    row = row.replace(substitution[0],substitution[1])
                    #print(row)

            
            if target_patterns:
                for substitution in target_patterns:
                            #        |target pattern|    |replacement sub|  
                    row = re.sub(r'' + substitution[0] + '', substitution[1], row)
                    #print(row) 

            if open_patterns!=[]:
                for substitution_tuple in open_patterns:
                    try:
                        to_be_replaced_pat=re.compile(r''+ substitution_tuple[0] +'')
                        replace_pat=re.compile(r''+ substitution_tuple[1] +'')
        
                        to_be_replaced_value=to_be_replaced_pat.match(row).group(1)
                        replace_value=replace_pat.match(row).group(1)

                        row=re.sub(to_be_replaced_value, replace_value, row)
                        #print(row)
                    except:
                        rows_not_modif_lst.append(row)
                        rows_not_modif_flag=True
                        # print(row)
                        # print('+'*100)
                        # a=input()
            
            #print(row)
            #a=input()    
            list_rows[ind] = row
            #print('\nlist_rows\n',list_rows)

        if rows_not_modif_flag==True:
            print('§§§§§§§§§§§ Some rows are not affected by the substitution §§§§§§§§§§§')

        if (rows_not_modif==False): return list_rows
        else: return list_rows, rows_not_modif_lst
    
    
    def FASTA_cpt_seq(self, list_rows):
        """
        Version: 1.0

        Name History: FASTA_cpt_seq

        This function compact the sequence after each FASTA header.
        It is possible that FASTA files downloaded from online sources could have
        the sequences splitted in multiple rows by '\n'. 
        The purpose is to arrange each sequence in a unique string.

        INPUT : list_rows   List[Str]   List of rows of a FASTA file. 
        OUTPUT: seq_compa   List[Str]   List of rows of a FASTA file with the sequences
                                        compacted in unique rows. 
        """    
        seq_compa=[]
        sequence_row=''

        for row in list_rows:    
            if row[0]=='>':                         
                if sequence_row!='':               # If the sequence line is empty but there is a FASTA header than this is the first header.
                    seq_compa.append(sequence_row) # Otherwise, it is a sequence that belongs to the current header and then it can be appended.
                seq_compa.append(row)              # If the current line is a FASTA header then append to the list.
                sequence_row=''                    # This means that you are expecting for a new sequence into the next rows.
            else:
                sequence_row += row  # Append the current part of sequence to the whole sequence.
        return  seq_compa


    def check_format(self, list_rows):
        """
        Version: 1.0

        Name History: check_format

        This function controls if the first row in a datafile contains the expected data.
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
            print('XXXXXX  WARNING - WRONG FILE FORMAT  XXXXXX') 

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
                row=row.replace('locus_tag='+target_tag, 
                                'gene:'+target_tag+' transcript:'+target_tag)
                                #gene:gene-   The first suggestion was to add gene:gene- but later was to remove them.

            FASTA_lst4PoGo.append(row)
                #print(row)
                #a=input()
        return FASTA_lst4PoGo    



    def protein_track(self, prot_list=[], strand='NC_006273.2', bed_fn='test1.bed'):
        """
        Version: 1.0

        Name History: protein_set_bed - protein_track

        Receive a list of protein codes and create the .bed track with the genomic locations of the protins.
        Example of multi exon protein in input:
          UL37 - [['52573', '53060', '-'], ['51302', '51344', '-'], ['50262', '51197', '-']]
        INPUT:  self.prot_pep_index
                self.prot_CDS_index
                self.prot_PSMint_index
        OUTPUT:
        """
        
        if prot_list==[]:
            prot_list=self.prot_pep_index.keys()

        print(f'Start processing {len(prot_list)} proteins')
        prot_track_fh = open(bed_fn,'w')

        if strand=='NC_006273.2': chr_name = 'NC_006273.2' # Find the chromosome name
        else:
            print('Option under construction. Check how to manage chromosome name in Human annotations\FASTA')

        BED_rows = []
        prot_row = ''

        for protein in prot_list:
            prot_row=''
            CDS_block=self.prot_CDS_index[protein] # Extract the block of CDS coordinates

            Strand=CDS_block[0][2]                 # The strand is in the third position of the first CDS features record

            if Strand =='+':                       # If strand positive
                chromStart=CDS_block[0][0]              # Chromosome Start is in the first position of the first record
                chromEnd=CDS_block[-1][1]               # Chromosome End is in the second position of the last record
            else:
                chromStart=CDS_block[0][1]              # In the negative case the bonduaries coordinates 
                chromEnd=CDS_block[-1][0]               # inverted their orders in the respective records

            print(self.prot_CDS_index[protein])
            print(f'{chromStart} - {chromEnd}')
            

            Score='1'
            tickStart = chromStart
            tickEnd = chromEnd
            blockCount = str(len(CDS_block))

            itemRGB = self.prot_PSMint_index[protein][2]  # Fetch the protein intensity into the proper index

            blockSizes_str = ''
            blockStarts_str = ''
            for CDS in CDS_block:
                print('CDSb',CDS_block)
                if Strand =='+':
                    CDS_start=int(CDS[0])
                    CDS_end=int(CDS[1])
                else:
                    CDS_start=int(CDS[1])
                    CDS_end=int(CDS[0])

                blockSizes_str+=str(abs(CDS_end-CDS_start)) + ','
                blockStarts_str+=str(abs(int(chromStart)-CDS_start)) + ','

            prot_row += (chr_name+'\t'+chromStart+'\t'+chromEnd+'\t'+protein+'\t'+Score+'\t'+
                         Strand+'\t'+tickStart+'\t'+tickEnd+'\t'+itemRGB+'\t'+blockCount+'\t'+
                         blockSizes_str+'\t'+blockStarts_str+'\n')

            print(prot_row)
            print('-'*100)


            BED_rows.append(prot_row)

        for row in BED_rows:
            prot_track_fh.write(row)
        prot_track_fh.close()

            #print(f'{protein} - {self.prot_index[protein]}')


    def PoGo_input_table(self, experiment_tag='exp', out_file_name=''):
        """
        Version: 1.0

        Name History: PoGo_input_table

        This function creates a peptide input table suitable for the PoGo software.
        Considering the format of the input table for Proteogenome, 
        this function removes the first column of the original input table. 
        Then, insert a new column into the first position with an experiment tag
        that will be reported on each peptide row.  

        INPUT : self.imput_table
        OUTPUT:
        """
        self.PoGo_it=np.delete(self.input_table,0,1)    # Remove the protein codes column. self.input_table[:,1::]
        self.PoGo_it=np.delete(self.PoGo_it,1,1)        # Remove the PTM column
        insert_experiment_name=np.full((len(self.PoGo_it),1), experiment_tag)  # Generate the experiment tag column.
        self.PoGo_it=np.concatenate((insert_experiment_name,self.PoGo_it), axis=1) # Insert the experiment tag column as a first column of the table.

        if out_file_name: 
            print('making PoGo input table')
            self.make_sep_file(out_file_name, self.PoGo_it, sep='\t') # Create the file with the the PoGo input table.


    def filter_peptides(self, PoGo_peptides):
        """
        Version: 1.0

        Name History: filter_peptides

        This function  filter the peptides mapped by PoGo. 
        The filter criteria is based on the proteins identified in the proteomics data 
        provided to the software.
        The peptides coordinates are comlpared with the CDS coordinates. 
        If the peptide coordinates between the CDS coordinates, then the peptides will 
        be included into the peptide map.  

        INPUT :
        OUTPUT:
        """
        row, col=PoGo_peptides.shape 

        # Initialise a smaller peptide table #
        coord_pep_strand=PoGo_peptides[:,1:6] # Slice the PoGo mapped peptide table considering only:
                                              # | start | end | peptide sequence | score | strand |
        coord_pep_strand=np.delete(coord_pep_strand,3,1)  # Delete the 'score' column

        print('Subsetting of PoGo_peptides\n',coord_pep_strand)

        protein_set=np.array([]) # Set of unique protein codes

        peptide_seq_array=coord_pep_strand[:,2]  # Consider all the peptide sequence
        for pep_seq in peptide_seq_array:
            protein_block= self.pep_prot_index[pep_seq] # Fetch the protein codes where the current peptide has been found
            protein_set=np.concatenate((protein_set,protein_block))
        protein_set=np.unique(protein_set)   # Shrink the protein code collection to the unique codes

        ################### GENERATE THE ALLOWED GENOMIC SPACE ###################
        allowed_genomic_space={}
        for unique_prot_code in protein_set:
            int_CDS_coordinates=[]                              
            for str_CDS_coord in self.prot_CDS_index[unique_prot_code]: # Considering the current protein code coordinates.
                from_str_to_int=list(map(int,str_CDS_coord[0:2]))       # Convert from Str to Int all the CDS coordinates.
                from_str_to_int.append(str_CDS_coord[-1])               # Append the strand tag in a Str format
                int_CDS_coordinates.append(from_str_to_int)             # Increase the coordinate bolck
            allowed_genomic_space[unique_prot_code]=int_CDS_coordinates # Append the new piece of allowed genomic coordinates.

        print(allowed_genomic_space)
        ##########################################################################

        for pep_row_index, peptide_row in enumerate(coord_pep_strand):    # Iterate over the PoGo peptide map
            peptide_coord_1=int(peptide_row[0])      # Fetch peptide genomic coordinates
            peptide_coord_2=int(peptide_row[1])
            peptide_sequence=peptide_row[2] 
            peptide_strand=peptide_row[3] 
            peptide_to_protein=self.pep_prot_index[peptide_sequence] # Fetch the set of proteins where the peptide has been found.

            # ---------- COORDINATES COMPARISON ---------- #  
            for protein in peptide_to_protein: # Iterate the set of protein 
                
                CDS_block=allowed_genomic_space[protein] # Fetch the genomic coordinates of the CDS of the protein where the peptide has been found.
                for CDS in CDS_block:
                    CDS_coord_1=CDS[0]
                    CDS_coord_2=CDS[1]
                    if peptide_strand=='+':
                        if (peptide_coord_1>=CDS_coord_1) and (peptide_coord_2<=CDS_coord_2): # VALID genomic coordinates
                            pass
                        else:                                                                 # INVALID genomic coordinates
                            PoGo_peptides=np.delete(PoGo_peptides,pep_row_index,0)     # REMOVE THE PEPTIDE FROM THE PoGo PEPTIDE TABLE
            
        print('PoGo_peptides')
        print(PoGo_peptides)
            # -------------------------------------------- #
        

# *************************************************************************************** #
# ********************************** Proteogenome DEMO ********************************** #

    def dummy_peptides(self, prot_sequence, pep_min_length, pep_max_length):
        """
        Version: 1.0

        Name History: dummy_peptides

        INPUT :

        OUTPUT:
        """

        prot_sequence = prot_sequence.replace('\n','')
        prot_sequence = prot_sequence.replace('K','|')  # LYSINE
        prot_sequence = prot_sequence.replace('R','|')  # ARGININE      
        
        new_peptides = prot_sequence.split('|') # The protein sequence is now a list of peptides
        
        peptides = []
        for pep in new_peptides:                # Filter the peptides with the lenght raange 
            if (len(pep) >= pep_min_length) & (len(pep) <= pep_max_length): 
                peptides.append(pep)   
        
        return peptides


    # def dummy_peptides_index(self, prot_sequences_FASTA_fn='', pep_min_length=4, pep_max_length=30):
    #     """
    #     Version: 1.0

    #     Name History: generate_peptides - dummy_peptides_index

    #     This function performs the synthetic trypsinisation of a protein sequence FASTA file.
    #     The output is a dictionary with the protein ID as key and the protein sequence as value.

    #     INPUT : prot_sequences_FASTA_fn   Str   File name of protein FASTA sequence file.
    #             pep_max_length  Int       Max lenght of the generated peptides.
    #             pep_min_length  Int       Min lenght of the generated peptides.
    #     OUTPUT: peptides        Dict      Key: [Str]    Protein ID.
    #                                       Val: [Str]    Protein sequence.       
    #     """

    #     #prot_ID_pat = re.compile(r'.*?protein=.*?\s(.*?)\sprotein_id=') # Pattern for protein ID
    #     prot_ID_pat = re.compile(r'.*?gene=.*?(.*?)\slocus_tag=') # Pattern for protein ID

    #     #FASTA_in_lst = self.file_to_lst(prot_sequences_FASTA_fn) # Upload the FASTA file with protein sequences in a list 
    #     FASTA_in_lst = self.FASTA_cpt_seq(self.FASTA_lst)         # Compact possible multilines sequences       

    #     self.dummy_input_matrix = np.array([[]])
    #     table_row=[]
    #     current_protein=True # This boolean signals if the data belong to the current protein. 
    #                          # This is used for combine the protein ID from the FASTA header with his sequence.

    #     for ind,prot_sequence in enumerate(FASTA_in_lst): # Iterate over the list that contains the FASTA file rows.
    #         if prot_sequence[0] == '>':# This condition is for the FASTA headers.
    #             # print(f'search header in row {ind}')
    #             # try:
    #             protein_ID=prot_ID_pat.match(prot_sequence).group(1) # Extract the PROTEIN ID
    #             #     print(protein_ID)
    #             # except:
    #             #     print('protein ID NOT FOUND --------------------')
    #             #     print(prot_sequence)
    #             #     a=input()
    #             current_protein=True     # Signals that the next rows belongs to this protein.
    #         else:                   # Generate the group of PEPTIDES that belong to the current protein.
    #             current_prot_pep=self.dummy_peptides(prot_sequence,pep_min_length, pep_max_length) 
    #             #print(f'the list of peptides for {protein_ID} is {current_prot_pep}')
    #             current_protein=False    # Signals that the next row will belong to a new protein.

    #         if not current_protein:
    #             # print(ind)
    #             # print(protein_ID)
    #             current_prot_pep.insert(0,protein_ID) # Put the protein ID before the peptide sequencies forming a new table row.
    #             # print(f'the complete row is {current_prot_pep}')
    #             # print(type(current_prot_pep))
    #             current_prot_pep=np.array(current_prot_pep,dtype=object) # Convert the list in a np.array.

    #             # Check if the array to concatenate has the same bumber of columns of the peptide matrix
    #             # pep_mat_rows=np.shape(self.dummy_input_matrix)[0] 
    #             pep_mat_cols=np.shape(self.dummy_input_matrix)[1]    # Number of columns in the protein matrix
    #             current_prot_pep_cols=np.shape(current_prot_pep)[1]  # Number of columns in the current Protein
    #             col_dif=pep_mat_cols-current_prot_pep_cols 

    #             if col_dif<0: # The number of peptides for this protein is greater than the number of peptides that belong to each of the proteins registered in the table so far


    #                 self.dummy_input_matrix=np.concatenate((self.dummy_input_matrix, [current_prot_pep]))
    #   # #         print(prot_sequence)
    #   #           prot_sequence = prot_sequence.replace('\n','')
    #   #           prot_sequence = prot_sequence.replace('K','|')  # LYSINE
    #   #           prot_sequence = prot_sequence.replace('R','|')  # ARGININE      
    #   #   #         print(prot_sequence)
                
    #   #           new_peptides = prot_sequence.split('|') # The protein sequence is now a list of peptides
    #   #   #         print(new_peptides)
                
    #   #           append_peptides = []
    #   #           for pep in new_peptides:
    #   #               if (len(pep) >= pep_min_length) & (len(pep) <= pep_max_length): # Filter the peptides with the lenght criteria
    #   #                   append_peptides.append(pep)
    #   #   #         print(append_peptides)    
                
    #   #           if len(append_peptides) > 0:
    #   #               peptides += append_peptides            
            
    # #         a=input()


    def dummy_input(self, peptides='', out_file_name='', dummy_exp_name = 'exp1',
                    prot_sequences_FASTA_fn='', pep_min_length=4, pep_max_length=30,
                    PSMs_range=1, pep_int_range=[100,10000]):
        """
        Version: 1.0

        Name History: dummy_PoGo_peptides - dummy_input   

        This function generates a file .txt that contains dummy peptides in a format suitable for PoGo software.

        INPUT : peptides            List[Str]       List of peptide sequence.
                out_file_name       Str             Name of the file that will contains the dummy peptides.
                dummy_exp_name      Str             Dummy experiment name.
        OUTPUT: dummy_input_matrix  np.array[Str]   The matrix with the dummy peptides.
        """

        prot_ID_pat = re.compile(r'.*?gene=.*?(.*?)\slocus_tag=') # Pattern for protein ID

        #FASTA_in_lst = self.file_to_lst(prot_sequences_FASTA_fn) # Upload the FASTA file with protein sequences in a list 
        FASTA_in_lst = self.FASTA_cpt_seq(self.FASTA_lst)         # Compact possible multilines sequences       

        dummy_input_matrix = np.array([['','','','','']])
        table_row=[]
        current_protein=True # This boolean signals if the data belong to the current protein. 
                             # This is used for combine the protein ID from the FASTA header with his sequence.

        for ind,prot_sequence in enumerate(FASTA_in_lst): # Iterate over the list that contains the FASTA file rows.
            #if ind>4: break # ################REMOVE THIS ROW IN THE FINAL VERSION################
            
            if prot_sequence[0] == '>':# If this row is a FASTA headers.
                protein_ID=prot_ID_pat.match(prot_sequence).group(1) # Extract the PROTEIN ID
                current_protein=True       # Signals that the next rows belongs to this protein.
            
            else:                   # Generate the group of PEPTIDES that belong to the current protein.
                current_prot_pep=self.dummy_peptides(prot_sequence,pep_min_length, pep_max_length) 
                current_protein=False      # Signals that the next row will belong to a new protein.

            if not current_protein:
                for new_peptide in current_prot_pep: # Iterate over the peptides generated from the current protein sequence.
                    
                    current_prot_pep=np.array([protein_ID,new_peptide],dtype=object) # Convert the list in a np.array.
                    

                # ***** RANDOM PSMs ***** #
                    if type(PSMs_range)==list:  
                                          # lower | upper
                                          # bound | bound
                        rand_PSM=randrange(PSMs[0],PSMs[1]) # Generate random PSMs
                    else:
                        rand_PSM=1
                    rand_PSM=np.array(['',rand_PSM])
                    #current_prot_pep=np.append(current_prot_pep,[rand_PSM],axis=0)
                    current_prot_pep=np.append(current_prot_pep,[rand_PSM],axis=0)

                # ***** RANDOM INTENSITY ***** #
                    rand_intensity=np.array(randrange(pep_int_range[0],pep_int_range[1]))
                    current_prot_pep=np.append(current_prot_pep,[rand_intensity],axis=0)
                    dummy_input_matrix=np.concatenate((dummy_input_matrix, [current_prot_pep]))
        
        dummy_input_matrix=np.delete(dummy_input_matrix,0,0) # Remove the first initialisation row
        

        #************* Output File creation
        if out_file_name !='':
            
            fh = open(out_file_name, 'w')
            for ind, row in enumerate(dummy_input_matrix):
                write_row = ''
        #         print(pep)
                write_row += dummy_exp_name + '\t' + row[1] + '\t' + str(row[2]) + '\t' + str(row[3]) + '\n'
                # write_row += dummy_exp_name + '\t' + pep + '\t' + str(PSM[ind]) + '\t' + str(inten[ind]) + '\n'
        #         print(row)
                fh.write(write_row)
        #         a=input()
            fh.close()

        return dummy_input_matrix
    
    def prot_not_represented(self, CDS_ind, pep_ind):
        """
        Version: 1.0

        History Name :

        This function cleans the self.prot_CDS_index class attribute of 
        protein codes that are not present in the self.prot_pep_index attribute class.
        One of the reason for the discrepancy could comes from the dummy_input function. 
        Here two thresholds have been set to define the length of the 
        dummy peptides (max - min).
        This could generate the case that some proteins are not represented in the 
        self.prot_pep_index attribute because none of its dummy peptides met the length 
        requirements set for the simulation..

        INPUT :
        OUTPUT:
        """
        for protein in list(CDS_ind.keys()): # Consider only the protein IDs
            if protein not in list(list(pep_ind.keys())):
                del CDS_ind[protein]

        # for protein in list(self.prot_CDS_index.keys()): # Consider only the protein IDs
        #     if protein not in list(list(self.prot_pep_index.keys())):
        #         del self.prot_CDS_index[protein]
        return CDS_ind


    def dummy_input_nparrays(self, peptides='', out_file_name='', dummy_exp_name = 'exp1',
                             prot_sequences_FASTA_fn='', pep_min_length=4, pep_max_length=30,
                             PSMs_range=1, pep_int_range=[100,10000]):
        """
        Version: 1.0

        Name History: dummy_PoGo_peptides - dummy_input   

        This function generates a file .txt that contains dummy peptides in a format suitable for PoGo software.

        INPUT : peptides        List[Str]   List of peptide sequence.
                out_file_name   Str         Name of the file that will contains the dummy peptides.
                dummy_exp_name  Str         Dummy experiment name.
        OUTPUT:
        """

        prot_ID_pat = re.compile(r'.*?gene=.*?(.*?)\slocus_tag=') # Pattern for protein ID

        #FASTA_in_lst = self.file_to_lst(prot_sequences_FASTA_fn) # Upload the FASTA file with protein sequences in a list 
        FASTA_in_lst = self.FASTA_cpt_seq(self.FASTA_lst)         # Compact possible multilines sequences       

        #self.dummy_input_matrix = np.array([['','','','']])
        self.dummy_input_matrix = np.empty((0,4)) 
        table_row=[]
        current_protein=True # This boolean signals if the data belong to the current protein. 
                             # This is used for combine the protein ID from the FASTA header with his sequence.

        for ind,prot_sequence in enumerate(FASTA_in_lst): # Iterate over the list that contains the FASTA file rows.
            #if ind>4: break # ################REMOVE THIS ROW IN THE FINAL VERSION################
            if prot_sequence[0] == '>':# This condition is for the FASTA headers.
                protein_ID=prot_ID_pat.match(prot_sequence).group(1) # Extract the PROTEIN ID
                current_protein=True       # Signals that the next rows belongs to this protein.
            else:                   # Generate the group of PEPTIDES that belong to the current protein.
                current_prot_pep=self.dummy_peptides(prot_sequence,pep_min_length, pep_max_length) 
                current_protein=False      # Signals that the next row will belong to a new protein.

            if not current_protein: # If the boolean for current protein is false means that is time to generate table rows.
                for new_peptide in current_prot_pep: # Iterate over the peptides generated from the current protein sequence.
                    
                    current_prot_pep=np.array([protein_ID,new_peptide],dtype=object) # Convert the list in a np.array.
                    
            # else:
            #     PSMs_col=np.ones(np.shape(self.dummy_input_matrix)[0]) # Otherwise the PSM value is 1 for all the peptides
                    
                    # print(current_prot_pep)
                    # print(self.dummy_input_matrix)
                    self.dummy_input_matrix=np.append(self.dummy_input_matrix, [current_prot_pep],axis=0)
        
        print('INPUT MATRIX')
        print(f'Type input matrix - {type(self.dummy_input_matrix)}')
        print(f'Shape input matrix - {np.shape(self.dummy_input_matrix)}')
        print(f'input matrix -\n{self.dummy_input_matrix}')  
        # ******* GENERATE PEPTIDE PSMs ******* #

        if type(PSMs_range)==list:   # lower | upper  | Generates an array of random values with the      
                                     # bound | bound  | same number of rows present in the input table.
            PSMs_col=np.random.randint(PSMs[0],PSMs[1],size=np.shape(self.dummy_input_matrix)[0])
        else:
            PSMs_col=np.ones(np.shape(self.dummy_input_matrix)[0]) # Otherwise the PSM value is 1 for all the peptides
        PSMs_col=np.reshape(PSMs_col,(-1,1)) # The unspecified value is inferred to be the number of rows
        # PSMs_col=np.array([PSMs_col])
        # PSMs_col=PSMs_col.T
        print('\n\n---------------------------------------------------------------')
        print('NEW COLUMN')
        print(f'Type new column - {type(PSMs_col)}')
        print(f'Shape new column - {np.shape(PSMs_col)}')
        print(f'new column - \n{PSMs_col}')
        self.dummy_input_matrix=np.concatenate((self.dummy_input_matrix, PSMs_col),axis=1)
        #self.dummy_input_matrix=np.append(self.dummy_input_matrix, PSMs_col,axis=1)


        # ******* GENERATE PEPTIDE INTENSTIES ******* #

        pep_inten_col=np.random.randint(pep_int_range[0],pep_int_range[1],
                                   np.shape(self.dummy_input_matrix)[0])

        self.dummy_input_matrix=np.concatenate((self.dummy_input_matrix, pep_inten_col), axis=1)

#********************** File creation
    #     PSM   = np.ones(len(peptides,), dtype=int)
    #     inten = np.random.randint(1000, 10000, len(peptides), dtype=int)
        
    #     fh = open(out_file_name, 'w')
    #     for ind, pep in enumerate(peptides):
    #         row = ''
    # #         print(pep)
    #         row += dummy_exp_name + '\t' + pep + '\t' + str(PSM[ind]) + '\t' + str(inten[ind]) + '\n'
    # #         print(row)
    #         fh.write(row)
    # #         a=input()
    #     fh.close()
    

    def Proteogenome_Demo(self, prot_quantity=10, pep_length_range=[], 
                          pep_int_range=[50,5000], PTM_ptg=0):
        """
        This function generates dummy peptide table
        INPUT :
        OUTPUT:
        """

    def make_sep_file(self, out_file_name, input_array, sep=''):
        """
        Version: 1.0

        Name History: make_tab_file - make_sep_file

        This function creates a file from an np.array/List with separated fields.
     
        Example: [[...],
                  [...],
                    .
                  [...]]

        The function retrieves the number of columns by looking at the number of columns stored 
        in the first row of the input list/array.
        The column's content will be written into the file separated by the value of the sep variable.

        INPUT : out_file_name   [Str]       The file path where the new file will be created
                input_array     [np.array]
                                    or      Array or Listof data to write in the output file
                                  [List]
                                        
        
        OUTPUT: The file with the input array content
        """
        
        out_file_hand = open(out_file_name, 'w')
        if type(input_array)==list and type(input_array[0])==list: # In this case is a list of lists
            number_of_columns=len(input_array[0])
        else: number_of_columns=1                                 # In this case is a simple list
        
        if 'numpy.ndarray' in str(type(input_array)): number_of_columns = input_array.shape[1]

        for row in input_array:
            out_row =''
            for col_ind, col_data in enumerate(row):
                out_row += col_data
                if (col_ind < number_of_columns): out_row += sep#'\t'
            out_row += '\n'
            out_file_hand.write(out_row)
        out_file_hand.close()
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
    # #HCMV_instance.protein_index(HCMV_instance.GFF3_to_df)
    # print(HCMV_instance.prot_index)
    # print(f'**done** \n')S

    print('Create .bed track for the whole protein set of the HCMV')

    # protein_list_test=['RL1', 'UL29', 'UL150A']
    # HCMV_instance.protein_track()
    # print('hello world ....... at the end!!!!!')

