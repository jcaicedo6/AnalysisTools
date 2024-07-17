import ROOT
import pandas as pd
import numpy as np
import concurrent.futures
from tqdm import tqdm

def getData(key:str, inputRootFiles:str, columns:list=None, cut=None, fraction=None, chunk_size=100000, max_workers=None):
    """
    :param key: Tree name of RootFiles.
    :param inputRootFiles: {path}/files.root
    :param columns: Variables of interest.
    :param cut: root-like cut to be applied.
    :param fraction: fraction of sampling.
    :return: Pandas dataframe.
    """
    def process_chunk(chunk_range):
        start_idx, end_idx = chunk_range
        pbar = tqdm(total=end_idx - start_idx, desc="Processing chunk", unit=" rows")
        
        if columns:
            chunk_df = df.Range(start_idx, end_idx).AsNumpy(columns=columns)
        else:
            chunk_df = df.Range(start_idx, end_idx).AsNumpy()
            
        if chunk_df:
            chunk_data = pd.DataFrame(chunk_df)
            pbar.update(len(chunk_data))  # Update progress bar with number of rows processed
            pbar.close()  # Close progress bar for this chunk
            return chunk_data
        pbar.close()  # Close progress bar for this chunk if no data
        return pd.DataFrame()
    
    # Create the RDataFrame with or without a cut
    if cut:
        df = ROOT.RDataFrame(key, inputRootFiles).Filter(cut)
    else: df = ROOT.RDataFrame(key, inputRootFiles)

    # Cache the dataframe if columns are specified
    if columns:
        df = df.Cache(columns) 

    # Get the number of entries in the dataframe
    num_entries = df.Count().GetValue()
    if num_entries == 0:
        print(f'Number of entries in {inputRootFiles} = 0')
        if columns: df = df.AsNumpy(columns=columns)
        else: df = df.AsNumpy()
        data = pd.DataFrame(df)
    
    # Split the dataframe into chunks
    else: 
        num_chunks = (num_entries + chunk_size - 1) // chunk_size  # Ceiling division
        chunk_ranges = [(i * chunk_size, min((i + 1) * chunk_size, num_entries)) for i in range(num_chunks)]

        # Use ThreadPoolExecutor to process chunks in parallel
        data_list = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            for chunk_range in chunk_ranges:
                print(f"Processing chunk {chunk_range[0]}-{chunk_range[1]}")
                chunk_data = process_chunk(chunk_range)
                if not chunk_data.empty:
                    data_list.append(chunk_data)


        if not data_list:
            raise ValueError("No data was collected in the chunks.")

        data = pd.concat(data_list, ignore_index=True)
    if fraction: 
        if fraction > 1:
            data = data.sample(frac=fraction, replace=True, random_state=1)
            data.reset_index(drop=True, inplace=True)
        else:
            data = data.sample(frac=fraction, random_state=1, ignore_index=True)
            data.reset_index(drop=True, inplace=True)

    return data
#---------#
# import DataProcessing as dp
# df = dp.getData(key='treename', inputRootFiles='../ntuple.root', columns=['Var_1', 'Var_2',...], cut='Var_1 > 5.2 && Var_2 > 0.2')

# Function to save data to csv
def save_to_csv(data, path:str, key:str=None):
    if key:
        print(f'Saving {key} to {path}')
        data.to_csv(path+f'{key}.csv', index=False)
    else:
        print(f'Saving to {path}')
        data.to_csv(path+'.csv', index=False)


# Function to calculate dynamic limits for histograms
def calculate_limits(data, margin=0.1):
    min_val = np.min(data)
    max_val = np.max(data)
    range_val = max_val - min_val
    return min_val - margin * range_val, max_val + margin * range_val

# Function to create histograms with dynamic limits
def create_hist(df, leafs, name, title, nbinsx):
    '''
    Create a ROOT histograms from Pandas DF
    
    :param df: Input Pandas DataFrame.
    :param leafs: List of column names of interest.
    :param name: Name of each varible.
    :param title: Title of each histogram.
    :param nbinsx: Number of bins.
    '''
    df_dic = {key: df[key].values for key in leafs}
    rdf = ROOT.RDF.MakeNumpyDataFrame(df_dic)
    
    hist = []
    for i, leaf in enumerate(leafs):
        # Calculate dynamic limits
        data = df[leaf].dropna()  # Drop missing values
        lower_limit, upper_limit = calculate_limits(data)
        
        # Create histograms
        hist.append(rdf.Histo1D((f"{name[i]}_{leafs[i]}", title[i], nbinsx, lower_limit, upper_limit), leaf))
    
    return hist

#--------#

#vhist = dp.create_hist(df=df, leafs=leafs, name=[x+"_sig" for x in name], title=leafs, nbinsx=nbins)

def BestCutsPandasDF(df, leafsToCut, leafs, cuts):
    """
    Filter the given DataFrame based on specified cuts and return a Pandas DataFrame.

    :param df: Input Pandas DataFrame.
    :param leafsToCut: List of column names to apply cuts.
    :param leafs: List of column names of interest.
    :param cuts: List of tuples specifying the cut ranges for each leafToCut.
    :return: Filtered Pandas DataFrame.
    """
    # Read the pandas dataframe in form of dictionary
    df_dic = {key: df[key].values for key in leafs}
    # write the dictionary with numpy arrays to a RDataframe
    rdf = ROOT.RDF.MakeNumpyDataFrame(df_dic)

    # Make the conditionals for the cuts
    conditions = []
    for leaf, cut_range in zip(leafsToCut, cuts):
        lower_cut = leaf + ">" + str(cut_range[0])
        upper_cut = leaf + "<" + str(cut_range[1])
        conditions.extend([lower_cut, upper_cut])

    # Apply cuts to the RDataFrame
    condition_str = " && ".join(conditions)
    rdf = rdf.Filter(condition_str)


    # Convert the filtered RDataFrame back to a Pandas DataFrame
    columns = rdf.AsNumpy(leafs)
    filtered_data = pd.DataFrame(columns)

    return filtered_data
#------#


# df_wcuts = BestCutsPandasDF(df=df, leafsToCut_lowm=leafsToCut_lowm, leafs=leafs, cuts=bestcuts)