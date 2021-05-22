#protein sequence:
protein_sequence = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"
 
#used for testing the example:
#protein_sequence = "WHGCITVYWMTV"
#making a dictionary for storing P(H),P(S) and P(T) values ((H and S are the alpha and beta propensities)
value_dict = {
"E": [1.53, 0.26, 0.44],  #Glu
"A": [1.45, 0.97, 0.57],  #Ala
"L": [1.34, 1.22, 0.53],  #Leu
"H": [1.24, 0.71, 0.69],  #His
"M": [1.20, 1.67, 0.67],  #Met
"Q": [1.17, 1.23, 0.56],  #Gln
"W": [1.14, 1.19, 1.11],  #Trp
"V": [1.14, 1.65, 0.30],  #Val
"F": [1.12, 1.28, 0.71],  #Phe
"K": [1.07, 1.28, 1.01],  #Lys
"I": [1.00, 1.60, 0.58],  #Ile
"D": [0.98, 0.80, 1.26],  #Asp
"T": [0.82, 1.20, 1.00],  #Thr
"S": [0.79, 0.72, 1.56],  #Ser
"R": [0.79, 0.90, 1.00],  #Arg
"C": [0.77, 1.30, 1.17],  #Cys
"N": [0.73, 0.65, 1.68],  #Asn
"Y": [0.61, 1.29, 1.25],  #Tyr
"P": [0.59, 0.62, 1.54],  #Pro
"G": [0.53, 0.81, 1.68]   #Gly
}
#contains values for checking each 6 pair for predicting helices
big_array1 = []
check_array1 = []
#contains values for checking each 5 pair for predicting strands
big_array2 = []
check_array2 = []
#print(value_dict["Glu"])
 
#these 2 functions are used for checking if the particular strings expand for a longer helix or strand:
 
def left_check(string_to_check, index_here, helix_or_strand):
  
   if index_here <= 0:
       return 0
   else:
       if value_dict[string_to_check[index_here]][helix_or_strand - 1] + value_dict[string_to_check[index_here+1]][helix_or_strand] + value_dict[string_to_check[index_here + 2]][helix_or_strand] + value_dict[string_to_check[index_here ]][helix_or_strand]>= 4:
           return left_check(string_to_check, index_here - 1, helix_or_strand)
       else:
           return index_here
 
def right_check(string_to_check, index_here, helix_or_strand):
  
   if index_here >= len(string_to_check) - 1:
       return len(string_to_check) - 1
 
   else:
       if value_dict[string_to_check[index_here]][helix_or_strand ] + value_dict[string_to_check[index_here - 1]][helix_or_strand] + value_dict[string_to_check[index_here - 2]][helix_or_strand] + value_dict[string_to_check[index_here + 1]][helix_or_strand] >= 4:
           return right_check(string_to_check, index_here + 1, helix_or_strand)
       else:
           return index_here
 
#making a new array for >1(helices)
for index in range(len(protein_sequence) - 5):
   big_array1.append(protein_sequence[index: index + 6])
 
for val in range(len(big_array1)):
   #print("break" + str(val))
   count_greater_than_one = 0
   for i in range(7):
       #print("test" + str(i))
       if count_greater_than_one >= 4:
           check_array1.append(big_array1[val])
           break
       else:
           if i < 6:
               if value_dict[big_array1[val][i]][0] >= 1:
              
                   count_greater_than_one += 1
                   #print(count_greater_than_one)
 
#similar process for strands          
for index in range(len(protein_sequence) - 4):
   big_array2.append(protein_sequence[index: index + 5])
 
for val1 in range(len(big_array2)):
   count_greater_than_one = 0
   for i in range(6):
       if count_greater_than_one >= 3:
           check_array2.append(big_array2[val1])
           break
       else:
           if i < 5:
               if value_dict[big_array2[val1][i]][1] >= 1:
                   count_greater_than_one += 1
 
#print(big_array2)
#print(value_dict[big_array2[0][0]])
 
#now to finallty predict if the the particular results in the arrays are helices and strands:
#here are the operations
helix_array = []
helix_indices = []
for trav in range(len(check_array1)):
   start_index = protein_sequence.find(check_array1[trav])
   #print("gg" + str(start_index))
   s_index = left_check(protein_sequence, start_index, 0)
   #print(s_index)
   e_index = right_check(protein_sequence, start_index + 5, 0)
   #print(e_index)
   helix_array.append(protein_sequence[s_index:e_index + 1])
   for i in range(s_index, e_index + 1):
       helix_indices.append(i)
 
strand_array = []
strand_indices = []
for trav in range(len(check_array2)):
   start_index = protein_sequence.find(check_array2[trav])
   s_index = left_check(protein_sequence, start_index, 1)
   e_index = right_check(protein_sequence, start_index + 4, 1)
   strand_array.append(protein_sequence[s_index:e_index + 1])
   for i in range(s_index, e_index + 1):
       strand_indices.append(i)
   #strand_array.append(protein_sequence[s_index:e_index + 1])
'''print(helix_array)
print(left_check("WHGCITVYWMTV",4, 0))
print(helix_indices)
print(check_array1)
print(big_array1)
#print(left_check("WHGCITVYWMTV", 0, 1))
print(helix_indices)
print(strand_indices)'''
 
helix_final = []
[helix_final.append(x) for x in helix_indices if x not in helix_final]
 
strand_final = []
[strand_final.append(x) for x in strand_indices if x not in strand_final]
 
#print(helix_final)
#print(strand_final)
 
#moving on to assignment of the final answer by making 2 seperate strings and then conflict resolving for the areas that overlap
 
helix_string = ""
strand_string = ""
for j in range(len(protein_sequence)):
   helix_string = helix_string + "_"
   strand_string = strand_string + "_"
 
for j in range(len(protein_sequence)):
   if j in helix_final:
       helix_string = helix_string[:j] + "H" + helix_string[j + 1:]
#print(helix_string)
 
for j in range(len(protein_sequence)):
   if j in strand_final:
       strand_string = strand_string[:j] + "S" + strand_string[j + 1:]
 
#print(strand_string)
 
#conflict resolving between_helix_string and strand_string:
final_answer = ""
i = 0
hsum = 0
ssum = 0
ctr=0
while(i < len(helix_string)):
   #print(i)
   #print(final_answer)
   if helix_string[i] == "_" and strand_string[i] == "_":
       final_answer = final_answer + "T"
       i += 1
   elif helix_string[i] == "H" and strand_string[i] == "_":
       final_answer = final_answer + "H"
       i += 1
   elif helix_string[i] == "_" and strand_string[i] == "S":
       final_answer = final_answer + "S"
       i += 1
   else:
       while((i  < len(helix_string) and (helix_string[i] == "H" and strand_string[i] == "S")) ):
           hsum += value_dict[protein_sequence[i]][0]
           ssum += value_dict[protein_sequence[i]][1]
           i += 1
           ctr+=1
       if hsum > ssum:
           final_answer=final_answer+("H"*ctr)
       else:
           final_answer=final_answer+("S"*ctr)
       hsum=0
       ssum=0
       ctr=0
print(protein_sequence)
print(final_answer)
part_1 = final_answer[0:50]
part_2 = final_answer[50:100]
part_3 = final_answer[100:150]
 
print(part_1)
print(part_2)
print(part_3)
