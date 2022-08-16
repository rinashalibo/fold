import numpy as np
import pandas as pd
import os
from glob import glob
from multiprocessing import cpu_count
import tempfile
import subprocess
from joblib import Parallel, delayed
from os.path import join as pjoin, basename, dirname
from tqdm import tqdm
import time

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    if type(seq) == str:
        reverse_complement = "".join(complement.get(base, base)
                                     for base in reversed(seq))
    elif type(seq) == list:
        reverse_complement = [complement.get(
            base, base) for base in reversed(seq)]
    elif type(seq) == np.ndarray:
        reverse_complement = np.array(
            [complement.get(base, base) for base in reversed(seq)])

    return reverse_complement


def get_pairing_score(f_ps, n_bases=10):
    (
        # read_name = f"
        #         {read_index}_
        #         {flow_in_read}_
        read_index,
        flow_in_read,
        seq_to_fold_len
    ) = basename(f_ps).split("_")[:-1]
    n = int(seq_to_fold_len)  # this is important and you should keep this in the metadata

    # read the results from the ps file
    lines = list()
    with open(f_ps) as fh:
        found_start = False
        for line in fh:
            line = line.strip()
            if not found_start:
                if "%start of base pair probability data" in line:
                    found_start = True
                continue
            if not line.endswith("ubox"):
                break
            lines.append(line.split()[:-1])

    # working out the results of the RNAfold to a local folding probability
    df = pd.DataFrame(lines, columns=["base", "base2", "p"]).astype(
        {"base": int, "base2": int, "p": float}
    )
    # print('RNAfold output:\n',df.head(),df.shape)
    df = pd.concat(
        (df[["base", "p"]], df[["base2", "p"]].rename(columns={"base2": "base"}))
    )
    # print('after concat:\n',df.head(),df.shape)

    df.loc[:, "p"] = df["p"] ** 2
    # print('after square:\n',df.head(),df.shape)
    df = df.groupby("base").sum().reset_index()
    # print('after group by and sum:\n',df.head(),df.shape)

    df.loc[:, "base"] = n + 1 - df["base"]
    # print(df.head(),df.shape,n)

    df = df.set_index("base").reindex(np.arange(1, 1 + n_bases)).fillna(0)
    # print(df.head(),df.shape)

    # adapt to your metadata
    df = (
        df.assign(
            read_index = read_index,
            flow_in_read = flow_in_read,
            seq_to_fold_len=int(seq_to_fold_len)
        )
            .reset_index()
            .set_index("read_index") #TODO: maybe change the indexing
            .sort_index()
    )

    # print(df.head())
    return df

def get_folding_results_1file(workdir, fasta_in):
    T = 37
    n_jobs = cpu_count() - 1
    fasta_out = fasta_in.replace(".fasta", ".rnafold.fasta")
    print(f"/home/ubuntu/ViennaRNA-2.5.1/src/bin/RNAfold -p -P DNA -T{T} -i{fasta_in} -o{fasta_out} --noconv --noPS -j{n_jobs}")
    rnafold_rc = subprocess.call(
        f"/home/ubuntu/ViennaRNA-2.5.1/src/bin/RNAfold -p -P DNA -T{T} -i{fasta_in} -o{fasta_out} --noconv --noPS -j{n_jobs}".split(),
        # f"/home/ubuntu/ViennaRNA-2.5.1/src/bin/RNAfold -p -P DNA -T{T} -i{fasta_in} --noconv --noPS -j{n_jobs}".split(),
        cwd=workdir
    )
    for f in glob(workdir + '*rnafold*'):
        os.remove(f)
    print('folding results in {}'.format(fasta_in))
    return

def get_folding_results(fasta_in, T, n_jobs):
    verbose = True
    progress_bar = False
    folding_results = {}
    with tempfile.TemporaryDirectory(prefix="/mnt/ramdisk/") as tmpdir:
        # fold
        # print('fold ...')
        fasta_out = fasta_in.replace(".fasta", ".rnafold.fasta")
        # tmpdir = '/data/work/220322/pycharm/smallfastas'
        rnafold_rc = subprocess.call(
            f"/usr/local/bin/RNAfold -p -P DNA -T{T} -i{fasta_in} -o{fasta_out} --noconv --noPS -j{n_jobs}".split(),
            # f"/usr/local/bin/RNAfold -p -P DNA -T{T} -i{fasta_in} -o{fasta_out} --noconv --noPS".split(),
            cwd=tmpdir,
        )
        # parse results
        # print('parse ...')
        psfiles = glob(pjoin(tmpdir, "*.ps"))[:]
        dfs = []
        dfs = pd.concat(
            Parallel(n_jobs=n_jobs)(
                delayed(get_pairing_score)(psfile)
                for psfile in psfiles
                #                     for psfile in tqdm(psfiles)
            )
        )
        # dfs_1 = dfs.query("base<=1")

    return dfs

def fold_and_save_single_file (workdir, fasta):

    return

def fold_and_save (workdir, fasta_merged):
    # fasta_merged = '/data/work/220322/pycharm/smallfastas/004777-X0024correction.fasta'
    # > split -d -l 200000 004777-X0024full.fasta merged.fasta.
    # fasta_merged = '/data/work/220322/pycharm/fastas/merged.fasta'

    # split the fasta_merged to chunks
    chunk_files = '{}.merged.fasta.'.format(fasta_merged.split('.fasta')[0])
    os.system('split -d -l 200000 {} {}'.format(fasta_merged,chunk_files))

    n_jobs_default = cpu_count() - 1
    T = 37
    folding_results_files = list()
    for fasta_chunk_file in tqdm(
        # list(glob(fasta_merged)[:]), desc=f"Iterating split fasta file, T={T}"
        list(glob(chunk_files + "*")[:]), desc = f"Iterating split fasta file, T={T}"
    ):
    # for fasta_chunk_file in list(glob(fasta_merged + ".*")[:2]):
    #     print(fasta_chunk_file)
    #     parquet_chunk_file = fasta_chunk_file.replace(".fasta","") + ".parquet"
    #     csv_chunk_file = fasta_chunk_file.replace(".fasta","") + ".csv"
        # print(parquet_chunk_file)
        # folding_results_files.append(parquet_chunk_file)
        # folding_results_chunk = get_folding_results(fasta_chunk_file,T,n_jobs_default)
        get_folding_results_1file(workdir, fasta_chunk_file)
        # folding_results_chunk.to_parquet(parquet_chunk_file)
        # folding_results_chunk.to_csv(csv_chunk_file)
        # print(parquet_chunk_file)
    return

def genome_key2seq( genomKey ):

    # This function translates a 1-D key series into a sequence of nucs
    seq_length = np.int64(np.sum(genomKey))
    seq_out = {}
    seq_out['seq'] = np.zeros( seq_length, dtype=np.int64 )
    seq_out['ind_orig'] = np.zeros( seq_length, dtype=np.int64 )
    seq_out['ind_last_hmer'] = np.zeros( seq_length, dtype=np.int64 )

    locs0 = np.where( genomKey>0 )[0]
    base_type = np.mod( np.arange(0,len(genomKey)), 4 )
    hmer_vec = np.int8(genomKey[locs0])
    cs = np.concatenate(([0,], np.cumsum(hmer_vec)), axis=0)
    max_hmer = np.max( hmer_vec )
    # go over all Hmer level, write all indexes together
    for i in range(1,max_hmer+1):
        locsh = np.where(hmer_vec==i)[0]
        locs_new = cs[locsh]  # starting location of this Hmer in the list of nuc's
        for j in range(0,i): # write info on all Hmer nuc's
            seq_out['seq'][ locs_new+j ] = base_type[locs0[locsh]]
            seq_out['ind_orig'][ locs_new+j ] = locs0[locsh]
            seq_out['ind_last_hmer'][locs_new+j] = locs_new+i-1

    # old, straight fwd code, slow
    # c = 0
    #seq_out1 = {}
    #seq_out1['seq'] = np.zeros( seq_length )
    #seq_out1['ind_orig'] = np.zeros( seq_length )
    #seq_out1['ind_last_hmer'] = np.zeros( seq_length )
    #for i in range(0,len(locs0)):
        #hmer = np.int8(genomKey[locs0[i]])
        #seq_out1['seq'][ c:c+hmer ] = base_type[locs0[i]]
        #seq_out1['ind_orig'][ c:c+hmer ] = locs0[i]
        #seq_out1['ind_last_hmer'][c:c+hmer] = c+hmer-1
        #c = c+hmer

    return seq_out


def get_key_mat():
    # loadDir = '/data/Runs/004777_1/s3files/'
    # offset_file = loadDir + 'RunID004777_1.Result.signal_file_offset.monitoring'
    # info_file = loadDir + 'RunID004777_1.monitor.info.csv'
    # reg_file = loadDir + 'RunID004777_1.Result.regr_sig.monitoring'
    # tk_file = loadDir + 'RunID004777_1.monitor.trueKey.key'

    loadDir = '/data/Runs/004777002_1/output/'
    offset_file = loadDir + 'RunID004777002_1.Result.signal_file_offset.monitoring'
    info_file = loadDir + 'monitor-004777002_1-gt.csv'
    reg_file = loadDir + 'RunID004777002_1.Result.regr_sig.monitoring'
    tk_file = loadDir + 'monitor-004777002_1-gt.key'

    print('load info...')
    info = pd.read_csv(info_file)
    print('load true key...')
    tk = np.fromfile(tk_file,dtype=np.int8).reshape(-1,416) # or 404
    print('load regressed...')
    regressed = np.fromfile(reg_file,dtype=np.int16)
    reg_key = np.round(regressed.reshape(-1, 440) / 100).astype(int) # or 428
    print('load offset...')
    offset = np.fromfile(offset_file, dtype=np.uint64)
    return reg_key, info, tk, offset, reg_file

def prepare_fasta_for_folding(key_mat,adapter_idx,fasta_name):
    # loop over reads in key mat
    #   loop over flows in read
    #       if flow > 0, apply key2seq + revcomp and write to fasta
    pb_adapter = "ATCACCGACTGCCCATAGAGAGCTGAGACTGCCAAGGCACACAGGGGATAGG"[5:]

    flow_array = np.array(['T', 'G', 'C', 'A'])
    # key_mat_aligned = key_mat[info['Bead Index'].values]
    key_mat_aligned = key_mat
    # random
    np.random.seed(0)
    # randomidx = np.random.choice(len(key_mat), 100, replace=False)
    randomidx = np.random.choice(len(key_mat), len(key_mat), replace=False)
    idx_name = fasta_name + ".idx"
    np.save(idx_name,randomidx)

    info_name = fasta_name+".info"
    n_reads_max = 100000000
    n_reads_now = 0
    with open(info_name, "w") as finfo:
        with open(fasta_name, "w") as ffasta:
            for read_index in randomidx:
                if (n_reads_now > n_reads_max):
                    break
                read = key_mat_aligned[read_index]
                adapter_index = adapter_idx[read_index]
                key_seq = read[:adapter_index]
                a = genome_key2seq(key_seq)
                ind_orig = np.unique(a['ind_orig'])
                ind_last_hmer = np.unique(a['ind_last_hmer'])
                seq = "".join(list(flow_array[a['seq']]))
                # direction = info['Direction'][read_index]

                idx_list = np.insert(len(seq) - 1 - ind_last_hmer, 0, len(seq))[:-1]
                idx_list_adapter = idx_list + len(pb_adapter)
                seq_final = pb_adapter + revcomp(seq)
                for j,i in enumerate(idx_list_adapter):
                    # print(read_index,len(seq_final[:i]),ind_orig[j],seq_final[:i])
                    read_name = f"{read_index}_{ind_orig[j]}_{len(seq_final[:i])}"
                    # write a fasta line: metadata and sequence
                    ffasta.write(f">{read_name}\n{seq_final[:i]}\n")
                    finfo.write(f"{read_index} {ind_orig[j]} {len(seq_final[:i])} {adapter_index}\n")
                n_reads_now += 1

    return

def read_folding_results_old (workdir):
    # list of parquet files
    print('read parquets...')
    folding_results_parquets = list(glob(pjoin(workdir,'*parquet*')))[:]
    # workdir
    print(folding_results_parquets)
    df_folding_results = pd.concat(
        (pd.read_parquet(f) for f in folding_results_parquets[:])
    ).sort_index()
    flows = []
    read_index = []
    folding_probability = [] # -1 * np.ones(len(df_base1[:].groupby(["cram_index", "flow"]).groups))
    # df_base1 = df_folding_results.query("base>1 & base<5")
    df_base1 = df_folding_results.query("base<4")
    print('loop...')
    for name, group in df_base1[:].groupby(["read_index", "flow_in_read"]):
        read_index.append(int(group.index[0]))
        flows.append(group["flow_in_read"].values[0])
        folding_probability.append(np.mean(group["p"].values))
    flows = np.array(flows)
    read_index = np.array(read_index)
    folding_probability = np.array(folding_probability)

    savedata = {}
    savedata['flows'] = flows
    savedata['read_index'] = read_index
    savedata['folding_probability'] = folding_probability
    results_file = pjoin(workdir,'savedata3flows')
    np.save(results_file,savedata)
    # plt.figure()
    # h1 = plt.hist(folding_probability[pattern_label == 'G'], bins='auto', histtype='step', density=True, log=True,label='G')
    # h2 = plt.hist(folding_probability[pattern_label == 'G1'], bins=h1[1], histtype='step', density=True, log=True,label='G->GG')
    # plt.show()
    #
    # read files

    # generate graph and save

    return results_file

def align_flow_order_with_first_flow_of_true_key(self):
    """Align the flow order with the first flow of the true key.
    The flow order parameter is defined with respect to the first flow of the preamble.
    i.e. ['TACG'] means the 'T' corresponds to the first preamble flow.
    Here we align the flow order with the first flow of the true key.
    i.e. ['TACG'] means the 'T' corresponds to the first flow of the true key
    return:
    aligned_flow_order - list of characters.
    """

    num_permutations = (self.params["NumPreamble"] + self.params["NumBarcode"] + self.params["NumUncertain"]) % len(self.params["Flow_order"])
    aligned_flow_order = self.params["Flow_order"].copy()
    for i in range(num_permutations):
        temp = aligned_flow_order.pop(0)
        aligned_flow_order.append(temp)

    return aligned_flow_order

def class_error(self, true_key, regressed_key, ind_filter, hmer=None):
    """Calculate the total class error, and the class error per base
    arguments:
    true_key - Numpy array. The true key of the beads we will use for our statistics.
    regressed_key - Numpy array. The regressed key of the beads we will use for our statistics
    ind filter - A boolean mask. Size of true key. True for flows we consider in the calculation.
    hmer - int. Perform calculation for a specific hmer.
    return:
    class_error_series - Pandas Series. Contains class error per base (and total).
    """
    class_error_dict = {}

    # Total number of class errors
    diff_mat = true_key - regressed_key

    # if ind_filter is None:
    #    ind_filter = np.ones(true_key.shape,dtype=np.bool)

    # Set boolean mask for a specific hmer
    if hmer is None:
        ind_hmer = np.ones(true_key.shape, dtype=bool)
    else:
        ind_hmer = true_key == hmer

    # Class error  - all bases combined
    indices = ind_filter & ind_hmer
    class_error_dict["Total (%)"] = 100 * np.mean(diff_mat[indices] != 0)  # Number of errors / total number of flows

    # Number of class errors per base
    aligned_flow_order = self.align_flow_order_with_first_flow_of_true_key()  # align the flow order with first flow of true key
    for nucl_base in aligned_flow_order:
        init_flow_index_for_base = aligned_flow_order.index(nucl_base)  # Index of nucl_base in the aligned flow_order
        diff_mat_for_nucl_base = diff_mat[:, init_flow_index_for_base::4]  # Diff mat for a specific nucl_base
        indices_for_nucl_base = (ind_filter[:, init_flow_index_for_base::4] & ind_hmer[:, init_flow_index_for_base::4])  # Indices for a specific nucl_base
        class_error_dict[f"{nucl_base} (%)"] = 100 * np.mean(diff_mat_for_nucl_base[indices_for_nucl_base] != 0)

    # Converting to Pandas
    # class_error_series = pd.Series(class_error_dict)
    # class_error_series.name = (str(hmer) if hmer is not None else 'Total')
    class_error_series = pd.DataFrame.from_dict(class_error_dict, orient="index").T
    class_error_series.index = ["Total"]
    return class_error_series


def build_folding_mat_from_npy (npy_file):

    data = np.ndarray.tolist(np.load(npy_file,allow_pickle=True))
    flows = data['flows'].astype(int)
    read_index = data['read_index']
    folding_probability = data['folding_probability']
    m = np.max(flows)
    n = np.max(read_index)
    folding_mat = np.zeros([n + 1, m + 1])
    folding_mat[read_index, flows] = folding_probability

    # info file
    offset_file = '/data/Runs/004777_1/folding_monitor_set/monitoring/RunID004777_1.Result.signal_file_offset.monitoring'
    info_file = '/data/Runs/004777_1/folding_monitor_set/csvs/RunID004777_1.monitor.info.csv'
    reg_file = '/data/Runs/004777_1/folding_monitor_set/RunID004777_1.Result.regr_sig.monitoring'
    tk_file = '/data/Runs/004777_1/folding_monitor_set/monitoring/RunID004777_1.monitor.trueKey.key'

    # match regressed and true key
    info = pd.read_csv(info_file)
    tk = np.fromfile(tk_file,dtype=np.int8).reshape(-1,404)
    regressed = np.fromfile(reg_file,dtype=np.int16).reshape(-1,428)
    offset = np.fromfile(offset_file, dtype=np.uint64)

    tkidx = np.searchsorted(offset, info.loc[:, "Bead Index"])
    reg_match_tk = np.round(regressed[tkidx]/100)
    folding_match_tk = folding_mat[tkidx[tkidx < 250000]]
    validtkmask = np.where(tk > -1)
    diffmat = np.abs(tk - reg_match_tk[:, 24:])
    ber = np.sum(diffmat[validtkmask])/len(validtkmask[0])

    print(data)

    return

def select_reads_and_get_adapter_position(key_mat, info, tk, offset,filter_file):
    f = {}
    np.random.seed(0)
    reads_in_tk = np.random.choice(len(tk),3)

    return

def filter_reads(regressed, info, tk, offset,filter_file):
    f = {}
    reads_in_tk = np.where((np.sum(tk == -65, axis=1) > 3) & (np.sum((tk > -1) & (tk < 13), axis=1) > 100))[0]
    # reads_in_tk = np.where((np.sum(tk == -65, axis=1) == 0) & (np.sum((tk > -1) & (tk < 13), axis=1) > 100))[0]
    reads_in_regressed = np.searchsorted(offset, info.loc[:, "Bead Index"].values[reads_in_tk])

    reg_match = regressed[reads_in_regressed][:,28:]
    tk_match = tk[reads_in_tk][:,4:]

    err_mat = (reg_match != tk_match) & ((tk_match > -1) & (tk_match < 13))
    idx_2errs = np.where(np.sum(err_mat, axis=1) < 2)[0]

    f['adapter_in_tk'] = reads_in_tk
    f['adapter_in_regressed'] = reads_in_regressed
    f['idx_2errs'] = idx_2errs
    np.save(filter_file,f)

    lines, rows = np.where(tk[reads_in_tk][idx_2errs] == -65)
    first_adapter_in_read = rows[np.unique(lines, return_index=True)[1]] + 24 # add 24 for tk and regressed difference in shape
    # first_adapter_in_read = (np.ones(len(tk[reads_in_tk][idx_2errs]))*np.shape(regressed)[1]).astype(int)

    return regressed[reads_in_regressed][idx_2errs],first_adapter_in_read

def read_folding_results(workdir):

    reg, info, tk, offset = get_key_mat()

    loaddir = workdir
    data = pd.read_csv(loaddir + 'RunID004777_1.Result.regr_sig.monitoring.fasta.output')
    print('get bead index list from output')
    randomidx = np.sort(np.load(loaddir + 'RunID004777_1.Result.regr_sig.monitoring.fasta.idx.npy'))

    z = []
    for i in data.values[:, 0]:
        z.append(i.split(' '))
    fold_p = np.array(z).transpose()[2,:].astype(float)

    tmp = []
    for i in np.array(z).transpose()[0, :]:
        tmp.append(i.split('_'))

    fold_flow = np.array(tmp).astype(int)[:,1]
    fold_bead_index = np.array(tmp).astype(int)[:,0]
    # sort
    sortidx = np.argsort(fold_bead_index)

    fold_flow = fold_flow[sortidx]
    fold_p = fold_p[sortidx]
    fold_bead_index = fold_bead_index[sortidx]

    folding_mat = np.zeros((len(randomidx),np.shape(reg)[1]))
    folding_mat_rows = np.cumsum(np.clip(np.diff(fold_bead_index,prepend=fold_bead_index[0]), 0, 1))
    folding_mat[folding_mat_rows, fold_flow] = fold_p

    f = np.load(loaddir + 'RunID004777_1.Result.regr_sig.monitoring.filter.npy',allow_pickle=True)
    adapter_in_tk = np.ndarray.tolist(f)['adapter_in_tk']
    adapter_in_regressed = np.ndarray.tolist(f)['adapter_in_regressed']
    idx_2errs = np.ndarray.tolist(f)['idx_2errs']

    reg_match = reg[adapter_in_regressed][idx_2errs][randomidx][:,28:]
    tk_match = tk[adapter_in_tk][idx_2errs][randomidx][:,4:]
    fold_match = folding_mat[:,28:]

    return folding_mat


flow_order = 'TGCA'
# workdir = '/data/Runs/004777_1/output/folding/'
# reg_file = '/data/Runs/004777_1/output/RunID004777_1.Result_ForGT.regr_sig.monitoring'
# info_file = '/data/Runs/004777_1/output/ForGT-004777_1-gt.csv'

workdir = '/data/Runs/004777002_1/fold_16_8_full/'
os.makedirs(workdir, exist_ok=True)
# info = pd.read_csv(info_file)
# reg_file = '/data/Runs/004777_1/s3files/RunID004777_1.Result.regr_sig.monitoring'
# reg_file = '/data/Runs/004777_1/monitor_files_v_4_7/RunID004777_1.Result.regr_sig.monitoring'


print('get key mat ... ')
key_mat, info, tk, offset, reg_file = get_key_mat()
fasta_file = workdir + reg_file.split('/')[-1]+'.fasta'
filter_file = workdir + reg_file.split('/')[-1]+'.filter'
print('filter reads by alignment results')
# select_reads_and_get_adapter_position(key_mat, info, tk, offset,filter_file)
key_mat,adapter_idx = filter_reads(key_mat, info, tk, offset,filter_file)
print('prepare fasta for folding ...')
prepare_fasta_for_folding(key_mat,adapter_idx,fasta_file)
start = time.time()
print('fold and save ... ')
# fold_and_save(workdir, fasta_file)
get_folding_results_1file(workdir, fasta_file)
end = time.time()
print('folding time : ', end - start)

# print('read folding results')
# read_folding_results(workdir)