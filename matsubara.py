import numpy as np

num_ferm_freq = int(input("Enter the number of Fermionic Matsubara frequencies you want in your active window: ")) #number of fermionic matsubara frequencies
Temp = float(input("Enter the electronic temp you want (in eV) :"))

w_b = list(range(-2*(num_ferm_freq-1), 2*num_ferm_freq, 2)) #bosonic matsubara frequncies
w_f = list(range(-(num_ferm_freq-1), num_ferm_freq, 2)) #fermionic matsubara frequencies 
w_b_pairs = [] #fermionic matsubara frequency pairs (wp,wk) s.t. wp-wk = wb
w_b_pairs_number = [] #number of fermionic matsubara freq pairs for each wb
T = 0.3 #electronic temperature

# Create a dictionary mapping values to their index (1-based)
index_map = {value: idx + 1 for idx, value in enumerate(w_f)}

# Define the index function
def f(x):
    return index_map.get(x, None)

    
#Find fermionic frequency pairs for each bosonic frequency
for i in w_b:
    temp_list = []
    for j in w_f:
        for k in w_f:
            if j-k == i:
                temp_list.append( (f(j),f(k)) )
    w_b_pairs.append(temp_list)
    w_b_pairs_number.append(len(temp_list))


w_b_p = [] #bosonic matsubara frequncies_pos
w_b_pairs_p = [] #fermionic matsubara frequency pairs (wp,wk) s.t. wp-wk = wb_pos
w_b_pairs_number_p = [] #number of fermionic matsubara freq pairs for each wb_pos


for i in range(len(w_b)):
    if w_b[i] >= 0:
        w_b_p.append(w_b[i])
        w_b_pairs_p.append(w_b_pairs[i])
        w_b_pairs_number_p.append(w_b_pairs_number[i])

#w_b = 2PinT
w_b_pos_T = [i*np.pi*T for i in w_b_p]
counter = 0
with open('non_neg_bosonic_frequencies', 'w') as BF:
    for w in w_b_pos_T:
        BF.write(str(w)+'\n')
with open('fermi_freq_pairs_len', 'w') as FFL:
    for i in range(len(w_b_p)):
        FFL.write(str(w_b_pairs_number_p[i])+'\n')
with open('fermi_freq_pairs', 'w') as FF:
    for i in range(len(w_b_p)):
        for j in range(len(w_b_pairs_p[i])):
                counter += 1
                FF.write(str(w_b_pairs_p[i][j][0]) +'     '+ str(w_b_pairs_p[i][j][1])+'\n')
with open('matsubara_params', 'w') as MAT:
    MAT.write(str(num_ferm_freq)+'\n')
    MAT.write(str(counter)+'\n')
    MAT.write(str(Temp))
