import numpy as np
def bianma_GA(xx0 = [0.07061393 ,0.02359821 ,0.52098946],min_a = [0., -1., 0.5],min_b = [2., 2., 2.],size_pop = 10,precision = 1e-4):

    xx0_list = xx0
    n_dim = len(xx0_list)
    min_a_list = min_a
    min_b_list = min_b
    #
    def to_bin(value, num):#
        bin_chars = ""
        value = int(value)
        temp = value
        for i in range(num):
            bin_char = bin(temp % 2)[-1]
            temp = temp // 2
            bin_chars = bin_char + bin_chars
        return bin_chars.upper()
    def to_2b_arry(value, num): #
        char_2b = to_bin(value, num)
        return  np.array(list(char_2b),dtype=int)
    ###
    lb, ub = np.array(min_a_list) * np.ones(n_dim), np.array(min_b_list) * np.ones(n_dim)
    #
    import scipy.stats as stats
    ####
    XX0_arry = np.zeros(shape=(size_pop, n_dim))
    XX0_arry_2b = np.zeros(shape=(size_pop, n_dim))
    for i in range(n_dim):
        sigma = 1
        X_tenp = stats.truncnorm((min_a_list[i] - xx0_list[i]) / sigma, (min_b_list[i] - xx0_list[i]) / sigma, loc=xx0_list[i], scale=sigma)
        a_temp = X_tenp.rvs(size_pop)
        XX0_arry[:, i] = a_temp
    XX0_arry[0,:] = np.array(xx0_list)
    #
    for i in range(n_dim):
        XX0_arry_2b[:, i] = (XX0_arry[:, i] - min_a_list[i])/ (min_b_list[i] - min_a_list[i])
    ###
    lb, ub = np.array(min_a_list) * np.ones(n_dim), np.array(min_b_list) * np.ones(n_dim)
    Lind_raw = np.log2((ub - lb) / precision + 1)
    Lind = np.ceil(Lind_raw).astype(int)  #
    Lind_sum = 0
    for i in Lind:
        Lind_sum += i
    #
    def b2_arry_gray2rv_b(gray_code):
        #
        size_pop_, len_gray_code = gray_code.shape
        aa = gray_code
        aaa = np.zeros(shape=(size_pop_,len_gray_code))
        for j in range(len_gray_code):
            if j != 0:
                aaa[:, j] = aa[:, j] - aa[:, j - 1]
            else:
                aaa[:, j] = aa[:, j]
        return aaa % 2
    #
    Chrom0 = np.zeros(shape=(size_pop,Lind_sum),dtype=int)
    Chrom0_start = {}
    for i in range(size_pop):
        Chrom0_start[i] = []
    for i in range(n_dim):
        len_gray_code = Lind[i]
        mask = np.logspace(start=1, stop=len_gray_code, base=0.5, num=len_gray_code)
        dx = mask[-1] #
        ndx_arry = np.around(XX0_arry_2b[:, i] / dx)
        aaa_temp = np.zeros(shape=(size_pop,len_gray_code))
        for j,j_v in enumerate(ndx_arry):
            aaa_temp[j,:] = to_2b_arry(j_v,len_gray_code) #
        bbb = b2_arry_gray2rv_b(aaa_temp)
        for j in range(size_pop):
            Chrom0_start[j] += list(bbb[j])
    for j in range(size_pop):
        Chrom0[j] = np.array(Chrom0_start[j])
    #
    def gray2rv(gray_code):
        # Gray Code to real value: one piece of a whole chromosome
        # input is a 2-dimensional numpy array of 0 and 1.
        # output is a 1-dimensional numpy array which convert every row of input into a real number.
        _, len_gray_code = gray_code.shape
        b = gray_code.cumsum(axis=1) % 2
        mask = np.logspace(start=1, stop=len_gray_code, base=0.5, num=len_gray_code)
        return (b * mask).sum(axis=1) / mask.sum()
    return Chrom0









