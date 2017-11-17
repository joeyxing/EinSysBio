#!/home/joey/Apps/SysBio2/bin/python
''' Author: Weizhou Xing
    System biology module 3
    Part A
'''
import numpy as np

MATCH = 2
MISMATCH = -1
GAP = -4


class return_object(object):
    def __init__(self):
        pass


def myArgmax(l):
    '''
    find all the max index
    i.e. myArgmax([1,2,3,4,5,5,5,3,2,1]) = [4,5,6]
    l: iterable
    '''
    m = np.max(l)
    r = []
    for i in range(len(l)):
        if l[i] == m:
            r.append(i)
    return np.array(r)


def global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

    for i in range(m):
        A[i, 0] = i*GAP
    for j in range(n):
        A[0, j] = j*GAP

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = np.max([A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP])
            else:
                A[i+1, j+1] = np.max([A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP])

    r = return_object()
    r.alignment = "global"
    r.A = A
    r.s = s
    r.t = t
    r.score = A[-1,-1]
    r.stop_list = [(m-1, n-1)]

    return r


def semi_global_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
    
    r = return_object()
    r.alignment = "semi"
    r.A = A
    r.s = s
    r.t = t

    max_r = np.max(A[-1,1:]) # max of last row
    argmax_r = myArgmax(A[-1,1:])+1 # argmax of last row
    max_c = np.max(A[1:-1,-1]) # max of last column (exclude (m-1, n-1))
    argmax_c = myArgmax(A[1:-1,-1])+1 # argmax of last column (exclude (m-1, n-1))
    r.stop_list = []
    if max_r >= max_c:
        r.score = max_r
        for i in range(len(argmax_r)):
            r.stop_list.append((m-1,argmax_r[i]))
    if max_r <= max_c: 
        r.score = max_c
        for i in range(len(argmax_c)):
            r.stop_list.append((argmax_c[i], n-1))

    return r


def local_alignment(s="", t=""):
    m = len(s) + 1
    n = len(t) + 1
    A = np.zeros([m, n], dtype=np.int32)

    for i in range(m-1):
        for j in range(n-1):
            if s[i] == t[j]:
                A[i+1, j+1] = max(0, A[i,j]+MATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)
            else:
                A[i+1, j+1] = max(0, A[i,j]+MISMATCH, A[i,j+1]+GAP, A[i+1,j]+GAP)

    r = return_object()
    r.alignment = "local"
    r.A = A
    r.score = np.max(A[:,:])
    l = np.where(A == r.score)
    r.stop_list = []
    for i in range(len(l[0])):
        r.stop_list.append((l[0][i], l[1][i]))
    # r.stop_list.append(np.unravel_index(np.argmax(A), A.shape))    
    r.s = s
    r.t = t

    return r


def trace_back(A, i, j, s, t, ei=0, ej=0, alignment=None):
    '''
    A: score matrix
    i: current index of str1 (row)
    j: current index of str2 (column)
    s: sequence 1
    t: sequecne 2
    ei: row of end point
    ej: column of end point
    return: list of element [insertion pos of s, insertion pos of t, end(start) point]
    '''
    if alignment == "global":
        if i==ei and j==ej: # and!
            insertion_pos = [[list(),list(),(i, j)]]
            return insertion_pos
    elif alignment == "semi":
        if i==ei or j==ej: # or!
            insertion_pos = [[list(),list(), (i, j)]]
            return insertion_pos
    # local alignment
    elif alignment == "local":
        # pass
        if i==ei and j==ej:
            insertion_pos = [[list(),list(), (i, j)]]
            return insertion_pos
    else:
        print "`trace_back` Error: unknow alignment method."
        return

    if s[i-1] == t[j-1]:
        cost = MATCH
    else:
        cost = MISMATCH

    insertion_pos = list()

    b0 = A[i,j-1]+GAP == A[i,j]     # b[0] from ->
    b1 = A[i-1,j-1]+cost == A[i,j]  # b[1] from  \
    b2 = A[i-1,j]+GAP == A[i,j]     # b[2] from  |

    if b0: # if previous is from left
        insertion_pos0 = trace_back(A,i,j-1,s,t,alignment=alignment)
        for trace in insertion_pos0:
            trace[0].append(i)
        insertion_pos = insertion_pos + insertion_pos0

    if b1: # if previous is from upper-left
        insertion_pos1 = trace_back(A,i-1,j-1,s,t,alignment=alignment)
        insertion_pos = insertion_pos + insertion_pos1

    if b2: # if previous if from upper
        insertion_pos2 = trace_back(A,i-1,j,s,t,alignment=alignment)
        for trace in insertion_pos2:
            trace[1].append(j)
        insertion_pos = insertion_pos + insertion_pos2

    # stop recursive for local alignment
    if not b0 and not b1 and not b2 and A[i,j] == 0:
        insertion_pos = [[list(),list(), (i, j)]] 

    return insertion_pos


def print_sequences(s, t, trace_list, stopi, stopj, count=0, insertSpaces=False):
    """
    s: sequence 0
    t: sequence 1
    trace_list: insertion_pos list
    stopi: end (right) index of sequence 0
    stopj: end (right) index of sequence 1
    count: number of different start point
    insertSpaces: global=false, semi-global=true, local=true
    """
    for k in range(len(trace_list)):
        trace = trace_list[k]
        # deal with sequence 0 (s)
        x = list()
        if trace[0]:
            x.append(s[:trace[0][0]])
            for i in range(len(trace[0])-1):
                x.append(s[trace[0][i]:trace[0][i+1]])
            x.append(s[trace[0][-1]:])
        else:
            x.append(s)
        sequence0 = "-".join(x)
        # deal with squence 1 (t)
        x = list()
        if trace[1]:
            x.append(t[:trace[1][0]])
            for i in range(len(trace[1])-1):
                x.append(t[trace[1][i]:trace[1][i+1]])
            x.append(t[trace[1][-1]:])
        else:
            x.append(t)
        sequence1 = "-".join(x)

        if trace[2][0] < trace[2][1]:
            sequence0 = " "*(trace[2][1] - trace[2][0]) + sequence0
        else:
            sequence1 = " "*(trace[2][0] - trace[2][1]) + sequence1

        bar = ""
        bar_length = (len(sequence0) - len(s) + stopi if len(sequence0) < len(sequence1)
                                     else len(sequence1) - len(t) + stopj)
        for i in range(bar_length):
            if sequence0[i] == sequence1[i]:
                bar = bar + "|"
            elif sequence0[i] == "-" or sequence1[i] == "-":
                bar = bar + "-"
            elif (sequence0[i] == " "
                  or sequence1[i] == " " 
                  or i < max(trace[2][:])):
                  #or i > min()):
                bar = bar + " "
            else:
                bar = bar + "x"

        print "\nAlignment{}:".format(count+k)
        print sequence0+"\n"+ bar + "\n" + sequence1


def special_semi_global_alignment(s="", t=""):
    # maximum of 10 gaps for each sequence
    m = len(s) + 1
    n = len(t) + 1
    max_gap = 10
    A = np.zeros([m, n],dtype=np.int32)
    # m, n, [gap for s, gap for t]
    A_gap = np.zeros([m,n,2], dtype=np.int8)
    
    for i in range(m-1):
        for j in range(n-1):

            if s[i] == t[j]:
                # MATCH
                cost = MATCH
            else:
                cost = MISMATCH

            # A[i,j]+cost is the first so that when equal cost occurs
            # always choose match or mismatch
            if np.argmax([A[i,j]+cost, A[i+1,j]+GAP, A[i,j+1]+GAP]) == 0:
                A[i+1,j+1] = A[i,j] + cost
                A_gap[i+1,j+1,0] = A_gap[i,j,0]
                A_gap[i+1,j+1,1] = A_gap[i,j,1]

            #### gap for "s" 
            elif np.argmax([A[i,j]+cost, A[i+1,j]+GAP, A[i,j+1]+GAP]) == 1:
                #  if "gap for s" <= max_gap, then, increase "gap for s" --
                if A_gap[i+1,j,1] + 1 <= max_gap:
                    A[i+1,j+1] = A[i+1,j] + GAP
                    A_gap[i+1,j+1,0] = A_gap[i+1,j,0] + 1
                    A_gap[i+1,j+1,1] = A_gap[i+1,j,1]

                # if "gap for s" > max_gap and "t+gap more than match" and "gap for t" <= max_gap
                # increase "gap for t" |
                elif A[i,j]+cost < A[i,j+1]+GAP and A_gap[i,j+1,1] + 1 <= max_gap:
                    A[i+1,j+1] = A[i,j+1] + GAP
                    A_gap[i+1,j+1,0] = A_gap[i,j+1,0]
                    A_gap[i+1,j+1,1] = A_gap[i,j+1,1] + 1

                # only match or mismatch direction is ok \
                else:
                    A[i+1,j+1] = A[i,j] + cost
                    A_gap[i+1,j+1,0] = A_gap[i,j,0]
                    A_gap[i+1,j+1,1] = A_gap[i,j,1]

            #### gap for "t"
            elif np.argmax([A[i,j]+cost, A[i+1,j]+GAP, A[i,j+1]+GAP]) == 2:
                # if "gap for t" <= max_gap, then, increase "gap for t"
                if A_gap[i,j+1,1] + 1 <= max_gap:
                    A[i+1,j+1] = A[i,j+1] + GAP
                    A_gap[i+1,j+1,0] = A_gap[i,j+1,0]
                    A_gap[i+1,j+1,1] = A_gap[i,j+1,1] + 1

                # if "gap for t" > max_gap and "s+gap more than (mis)match" and "gap for s" <= max_gap
                # increase "gap for s"
                elif A[i,j]+cost < A[i+1,j]+GAP and A_gap[i+1,j,0] + 1 <= max_gap:
                    A[i+1,j+1] = A[i+1,j] + GAP
                    A_gap[i+1,j+1,0] = A_gap[i+1,j,0] + 1
                    A_gap[i+1,j+1,1] = A_gap[i+1,j,1]

                # only match or mismatch direction is ok
                else:
                    A[i+1,j+1] = A[i,j] + cost
                    A_gap[i+1,j+1,0] = A_gap[i,j,0]
                    A_gap[i+1,j+1,1] = A_gap[i,j,1]
    # print A[:,:]
    # print A_gap[:,:,0]
    # print A_gap[:,:,1]
    # r1 = [i for i in np.where(A_gap[-1,1:,0]<=max_gap)[0]
    #          if i in np.where(A_gap[-1,1:,1]<=max_gap)[0]]
    # r2 = [i for i in np.where(A_gap[1:-1,-1,0]<=max_gap)[0]
    #          if i in np.where(A_gap[1:-1,-1,1]<=max_gap)[0]]
    # print r1
    # print r2
    print np.where(A_gap[:,:,0] > 10)
    print np.where(A_gap[:,:,1] > 10)

    return max(max(A[-1,1:]), max(A[1:-1,-1]))


def read_fasta(path="ebolasequences-1.fasta"):
    with open(path) as f:
        fasta = f.read()

    sequence_list = fasta.split("\n>")
    for i in range(len(sequence_list)):
        start_index = sequence_list[i].find("\n")
        sequence_list[i] = sequence_list[i][start_index+1:]
        sequence_list[i] = sequence_list[i].replace("\n", "")

    return sequence_list


def align_trace_print(s, t, method=""):
    if method == "global":
        r = global_alignment(s=s, t=t)
    elif method == "semi":
        r = semi_global_alignment(s=s, t=t)
    elif method == "local":
        r = local_alignment(s=s, t=t)
        # TODO
        all_trace_list = list()
    else:
        print "Error! Unknown alignment method."
        return

    print "Score:", r.score

    n = 0
    for i in range(len(r.stop_list)):
        trace_list = trace_back(r.A,
                                r.stop_list[i][0], r.stop_list[i][1],
                                r.s, r.t,
                                alignment=r.alignment)
        if method == "local":
            # TODO only use the longest subsequence
            for x in trace_list:
                all_trace_list.append(x+
                    [(r.stop_list[i][0], r.stop_list[i][1])])
            continue
        print trace_list
        print_sequences(r.s, r.t,
                        trace_list, 
                        r.stop_list[i][0], r.stop_list[i][1],
                        count=n, insertSpaces=True)
        n = n + len(trace_list)
    
    if method == "local":
        # Only print the longest sequence(s)
        max_subseq_len = 0
        final_list = []
        for x in all_trace_list:
            subseq_len = max(x[3][1]-x[2][1],
                             x[3][0]-x[2][1])
            if subseq_len > max_subseq_len:
                max_subseq_len = subseq_len
                final_list = [x]
            elif subseq_len == max_subseq_len:
                final_list.append(x)
        for x in final_list:
            print_sequences(r.s, r.t,
                            [x],
                            x[3][0], x[3][1],
                            count=0)


def main():
    print "System Biology Module 3 Part A"
    
    # Part A.(a):
    print "\n(a): Global Alignment"
    str1 = "ACAAGGA"
    str2 = "ACAGG"
    align_trace_print(str1, str2, method="global")
    
    # Part A.(b):
    print "\n(b): Semi-global Alignment"
    str3 = "AGCCAATTACCAATTAAGG"
    str4 = "CCAATT"
    align_trace_print(str3, str4, method="semi")

    # Part A.(c):
    print "\n(c): Local Alignment"
    str5 = "AGCCTTCCTAGGG"
    str6 = "GCTTCGTTT"
    align_trace_print(str5, str6, method="local")

    # Part A.(d):
    print "(d):"
    [str7,str8] = read_fasta(path="C:/Users/xwz20/Desktop/EinSysBio/Module3/PartA/ebolasequences-1.fasta")
    print len(str7)
    print len(str8)
    print special_semi_global_alignment(s=str7,t=str8)


def test():
    # test
    str1 = "ACAAGGAATGTACATCATTAGCTAGTCA"
    str2 = "ACAGGACTATCGATGTCGATGCTAGT"
    align_trace_print(str1, str2, method="global")

    # test
    str3 = "AGCCAATTACCAATTAAGGCCTACATT"
    str4 = "CCAACACTT"
    align_trace_print(str3, str4, method="semi")

    # test
    str5 = "AGCCTTCTACGCTAGGG"
    str6 = "GCTTCGTACGTTT"
    align_trace_print(str5, str6, method="local")

    str7 = "AGCCAATTACCAATTAAGGCCTACATT"
    str8 = "CCAACACTT"
    print special_semi_global_alignment(s=str7,t=str8)


if __name__ == "__main__":
    # TODO: improve docstring
    main()

