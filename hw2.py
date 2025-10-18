def calculate_SoP(input_path, score_path, gopen, gextend):
    import numpy as np
    import pandas as pd 
    seq = pd.read_csv(input_path, header = None, comment = ">")
    score_mat = pd.read_csv(score_path, header = 0, comment = "#", sep = "\s+")
    #Part 1: Score Matrix 
    #Extend the score_mat to include gap openings(-) & extensions(*)
    score_mat.loc["*"] = gextend
    score_mat.loc[:,"*"] = gextend
    score_mat.loc["-"] = gopen
    score_mat.loc[:,"-"] = gopen
    #Used trial & error to find out that score(-,*) = gap extension penalty
    score_mat.loc["-","*"] = gextend
    score_mat.loc["*","-"] = gextend
    score_mat["a.a.1"] = score_mat.index
    score_mat = score_mat.melt(id_vars = "a.a.1", var_name = "a.a.2", value_name = "score")
    score_mat = score_mat.set_index(["a.a.1", "a.a.2"])
    #Part 2: Data Wrangling
    #Convert the fasta file into a pandas dataframe, with the rows as sequences & the columns as protein bases
    seq = seq.apply(list, axis = 1)
    seq = seq.apply(str)
    seq = seq.apply(list)
    seq = seq.apply(pd.DataFrame)
    seq = pd.concat(list(seq), axis=1).T
    seq = seq.iloc[:,2:(seq.shape[1]-2)]
    # n_seq: number of sequences
    # seq_len: length of the sequence alignment (including gaps/spaces) 
    [n_seq, seq_len] = seq.shape
    #Find the positions of the gaps
    gap = pd.DataFrame.to_numpy(seq == "-")
    #Gaps are either openings or extensions; gap_opening = gap - gap_extend
    gap_extend = gap * np.concatenate((np.zeros((n_seq, 1)), gap[:,0:seq_len - 1]), axis = 1) 
    #Mark the gap openers with dashes(-) and gap extenstions with stars(*)
    seq = seq.where(gap_extend == 0, other = "*").T
    #Part 3: Define a funtion to calculate the score between a pair of sequences
    def sop(i, j):
         seq_subset = seq.iloc[:,[i,j]]
         seq_subset.columns = (["a.a.1", "a.a.2"])
         seq_subset = seq_subset.set_index(["a.a.1", "a.a.2"])
         seq_subset = list(seq_subset.index)
         pair_score = score_mat.reindex(seq_subset)
         pair_score = sum(pair_score["score"])
         return pair_score
   #Part 4: Find all possible pairs of sequences & sum up their scores
    i = np.repeat(np.arange(n_seq), n_seq)
    j = np.arange(n_seq)
    j = j[:, np.newaxis].T
    j = np.repeat(j, n_seq, axis = 0)
    j = np.reshape(j, n_seq**2)

    msa = pd.DataFrame([i, j]).T
    msa = msa.rename(columns={0: "i", 1: "j"})
    msa = msa.loc[msa["i"] < msa["j"]]
    msa["score"] = msa.apply(lambda x: sop(x["i"], x["j"]), axis = 1)
    final_score = sum(msa["score"])
    return final_score




