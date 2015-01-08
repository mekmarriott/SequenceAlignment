#!/usr/bin/python
#!/usr/bin/env python

import sys
import sets

input_file = sys.argv[1]
output_file = sys.argv[2]
f = file(input_file)
out = open(output_file, 'w+')

#initializing all global variables
seqA = f.readline() #make in to string or array of char
seqB = f.readline() #make in to string or array of char
num_A = len(seqA)
num_B = len(seqB)
align_type = int(f.readline()) # 0 or 1
gap_penalty = f.readline().split() #array [gap_penalty_A, gap_ext_A, gap_penalty_B, gap_ext_B]
gap_penalty_A = float(gap_penalty[0])
gap_ext_A = float(gap_penalty[1])
gap_penalty_B = float(gap_penalty[2])
gap_ext_B = float(gap_penalty[3])
num_letters_seqA = int(f.readline()) 
letters_seqA = f.readline() 
num_letters_seqB = int(f.readline())
letters_seqB = f.readline() 

#create map for scores
score_map = {}
for line in f:
	arr = line.split()
	if len(arr) > 3:
		score_map[arr[3],arr[2]] = float(arr[4])

#print score_map
#initialize set for completed alignments
alignment_set = set()

#initialize alignment matrices (am, ix, iy)
am = [ [ [0.0,[False,False,False]] for j in range(num_A) ] for i in range(num_B) ] #alignment matrix index --> [seqB_index][seqA_index]
ix = [ [ [0.0,[False,False]] for j in range(num_A) ] for i in range(num_B) ] #Ix matrix index --> [seqB_index][seqA_index]
iy = [ [ [0.0,[False,False]] for j in range(num_A) ] for i in range(num_B) ] #Iy matrix index --> [seqB_index][seqA_index]

def populate():
	for i in range(1,num_B):
		for j in range(1,num_A):
			m_score = score_map[seqB[i-1], seqA[j-1]]
			m_score_match = am[i-1][j-1][0] + m_score
			m_score_x = ix[i-1][j-1][0] + m_score
			m_score_y = iy[i-1][j-1][0] + m_score

			m_max_score = round(max(m_score_match, m_score_x, m_score_y),5)

			if (round(m_max_score,5) > 0 or align_type == 0):
				am[i][j][0] = m_max_score
				if (i > 1 and j > 1):
					if round(m_max_score,5) == round(m_score_match,5):
						am[i][j][1][0] = True
					if round(m_max_score,5) == round(m_score_x,5):
						am[i][j][1][1] = True
					if round(m_max_score,5) == round(m_score_y,5):
						am[i][j][1][2] = True
				else:
					am[i][j][1] = [True, False, False]

			#for ix______________________________________________
			x_score_match = am[i][j-1][0] - gap_penalty_B
			x_score_x = ix[i][j-1][0] - gap_ext_B

			x_max_score = round(max(x_score_match, x_score_x),5)

			if (round(x_max_score,5) > 0 or align_type == 0):
				ix[i][j][0] = x_max_score
				if (j > 1):
					if round(x_max_score,5) == round(x_score_match,5):
							ix[i][j][1][0] = True
					if round(x_max_score,5) == round(x_score_x,5):
							ix[i][j][1][1] = True
				else:
					ix[i][j][1] = [False, True]

			#for iy______________________________________________
			y_score_match = am[i-1][j][0] - gap_penalty_A
			y_score_y = iy[i-1][j][0]- gap_ext_A

			y_max_score = round(max(y_score_match, y_score_y),5)

			if (round(y_max_score,5) > 0 or align_type == 0):
				iy[i][j][0] = y_max_score
				if (i > 1):
					if round(y_max_score,5) == round(y_score_match,5):
						iy[i][j][1][0] = True
					if round(y_max_score,5) == round(y_score_y,5):
						iy[i][j][1][1] = True
				else:
					iy[i][j][1] = [False, True]


def trace(matrix, i, j, alignmentA, alignmentB):
	if (i > 0 and j > 0):
		#print 'looking at ' + str(i) + ' ' + str(j) + ' in ' + str(matrix)
		#Traceback for match matrix______________________________________________________________________
		if matrix == 0:
			if (am[i][j][1][0]== True):
				trace(0,i-1,j-1, seqA[j-1] + alignmentA, seqB[i-1] + alignmentB) #making match
			if (am[i][j][1][1] == True):
				trace(1,i-1,j-1, seqA[j-1] + alignmentA, seqB[i-1] + alignmentB) #make shift in x, i.e. put gap in sequence B and shift one in sequence A
			if (am[i][j][1][2] == True):
				trace(2,i-1,j-1, seqA[j-1] + alignmentA, seqB[i-1] + alignmentB) #make shift in y, i.e. put gap in sequence A and shift one in sequence B
			if (am[i][j][1] == [False,False,False] and align_type == 1):
				alignment_set.add((alignmentA, alignmentB))
		elif matrix == 1:
			if (ix[i][j][1][0]== True):
				trace(0,i,j-1, seqA[j-1] + alignmentA, '_'+ alignmentB) #go back to match matrix, don't change anything
			if (ix[i][j][1][1] == True):
				trace(1,i,j-1, seqA[j-1] + alignmentA, '_'+ alignmentB) #extend gap in x, i.e. put one more gap in sequence B and shift one in sequence A
		elif matrix == 2:
			if (iy[i][j][1][0]== True):
				trace(0,i-1,j, '_' + alignmentA, seqB[i-1] + alignmentB) #go back to match matrix, don't change anything
			if (iy[i][j][1][1] == True):
				trace(2,i-1,j, '_' + alignmentA, seqB[i-1] + alignmentB) #extend gap in x, i.e. put one more gap in sequence A and shift one in sequence B
	else:
		alignment_set.add((alignmentA, alignmentB))

#SCRIPT_________________________________________________________

populate()
max_coordinates = []
max_score = 0.0

if align_type == 0:
# Calculate maximum point for GLOBAL MAX
	for i in range(1,num_B):
		if round(am[i][num_A-1][0],5) > round(max_score,5):
			max_score = am[i][num_A-1][0]
			max_coordinates = [[0, i, num_A-1]]
		elif round(am[i][num_A-1][0],5) == round(max_score,5):
			max_coordinates += [[0, i, num_A-1]]

	for j in range(1,num_A-1):
		if round(am[num_B-1][j][0],5) > round(max_score,5):
			max_score = am[num_B-1][j][0]
			max_coordinates = [[0, num_B-1, j]]
		elif round(am[num_B-1][j][0],5) == round(max_score,5):
			max_coordinates += [[0, num_B-1, j]]
else:
#Calculate maximum point for LOCAL MAX__________________________________
	for i in range(1,num_B):
		for j in range(1,num_A):
			if round(am[i][j][0],5) > round(max_score,5):
				max_score = am[i][j][0]
				max_coordinates = [[0,i,j]]
			elif round(am[i][j][0],5) == round(max_score,5):
				max_coordinates += [[0,i,j]]


#Print/Output RESULTS_____________________________________________________
print str(max_score)
out.write(str(max_score))
for points in max_coordinates:
	trace(points[0], points[1], points[2], "", "")
for pair in alignment_set:
	out.write('\n')
	out.write(pair[0])
	out.write(pair[1])
	print '\n'
	print pair[0]
	print pair[1]

		



