#Time Complexity: O(nm)
#Space Complexity: O(nm)
INDEL = 5
MATCH = -3
SUB = 1
def unrestricted(seq1, seq2, align_length):
    
    n = len(seq1)
    m = len(seq2)
    table = [[0 for i in range(n+1)] for j in range(m+1)] #Zero the table
    # back_pointers  = [None * align_length]
    # print(table)
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0: #top or left-most column
                table[i][j] = (INDEL * i) + (INDEL * j)
            else:
                matched = seq1[j-1] == seq2[i-1]
                #take the min score from the surrounding
                if matched:
                    min = table[i-1][j-1] + MATCH #diagnal cell - Last priority in case of tie
                else:
                    min = table[i-1][j-1] + SUB

                if min >= table[i-1][j] + INDEL: #top cell - Second priority in case of tie
                    min = table[i-1][j] + INDEL

                if min >= table[i][j-1] + INDEL: #left cell - Top priority in case of tie
                    min = table[i][j-1] + INDEL
                # print(min)
                table[i][j] = min
                
    print(table)
    return table[m][n] #final score

print(unrestricted("Hello", "Help", 5))
