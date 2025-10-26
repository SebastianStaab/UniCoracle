# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:58:03 2024

@author: JohnDoe2Go
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:35:53 2022

@author: JohnDoe
"""
from unicor import unicorp
import coracle
import pandas as pd
import time


def unicoracle(x, y, tax,
        n_features=100,
        uc=True,
        uc_threshold=None,          # XOR with uc_top_k
        uc_top_k=None,              # XOR with uc_threshold
        uc_method='pearson',          # 'pearson' or 'spearman'
        uc_transformation='relative_abundance',       # 'raw', 'relative_abundance', or 'clr'
        seed=None,        # <- drives determinism end-to-end
        coracle_kwargs=None):
    """
    UniCoracle is a fully automated analytical framework that combines UniCorP’s 
    bottom-up propagation approach 
    (https://academic.oup.com/ismecommun/advance-article/doi/10.1093/ismeco/ycaf174/8269840)
    with a subsequent and newly developed top-down skimming (TDS) approach employing the 
    Coracle ML framework 
    (https://academic.oup.com/bioinformatics/article/40/1/btad749/7484655) 
    to leverage inherent taxonomic structures of microbiome community data 
    (e.g., 16S rRNA gene amplicon sequencing) for increasing predictive accuracy 
    and reducing runtimes.
    
    Parameters
    ----------
    x : pd.dataframe
        ASV dataset
    y : pd.dataframe
        Continuous target variable
    tax : pd.dataframe
        Taxonomic hierarchy
    n_features : int, optional
        Number of features to select at each level. Default is 100.
    uc: bool, optional
        True if unicor should be applied, False if only the top down skimming approach should be used. 
    uc_threshold : float, optional
        UniCor metric threshold. Can only be between 0 and 1. Optimal values depend on the dataset. The default is 0.15.
    uc_top_k : int, optional
        Number of UniCorP selections per level (mutually exclusive with `uc_threshold`).
        Use when you want stable dimensionality and runtime. 
    uc_method : {'pearson', 'spearman'}, optional
        Correlation used in the UniCor metric. Default is 'pearson'.
        'spearman' is more robust to monotonic nonlinearity and outliers.
    uc_transformation : {'raw', 'relative_abundance', 'clr'}, optional
        Transformation for UniCor scoring only (the pipeline keeps counts for modeling).
        Default is 'relative_abundance'. 'clr' applies centered log-ratio.
    seed : int, optional
        Random seed propagated to Coracle (and underlying models) for repeatability.
        Set to an integer for deterministic runs; `None` uses estimator defaults.
    coracle_kwargs : dict, optional
        Extra options forwarded to Coracle during TDS and the final pass (e.g., models,
        CV scheme, scoring, n_jobs). See Coracle docs for accepted keys/values.
        
    Raises
    ------
    TypeError
        check types of input
    ValueError
        check allowed ranges for the thresholds, dimensions for the dataframes

    Returns
    -------
    resultO : pd.dataframe
        Final Coracle result on the lowest level
 """       

    
    
    ### check input
    ###########################################################################
    # check types
    if not isinstance(x, pd.DataFrame):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for ASV")
    if not isinstance(y, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe or series for target")
    if not isinstance(tax, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for tax")

    if not isinstance(n_features, int) or n_features <= 1:
        raise ValueError("ValueError exception thrown. Expected n_features to be a positive integer. In order to keep compositionality it has to be above 1.")

    if not isinstance(uc, bool):
        raise TypeError("TypeError exception thrown. Expected bool for uc")

    # uc_threshold / uc_top_k mutual exclusivity and typing
    if uc:
        if (uc_threshold is not None) and (uc_top_k is not None):
            raise ValueError("ValueError exception thrown. Use either uc_threshold or uc_top_k, not both.")

        if (uc_threshold is None) and (uc_top_k is None):
            # sensible default to keep dimensionality/runtime predictable
            uc_top_k = n_features

        if uc_threshold is not None:
            if not isinstance(uc_threshold, (int, float)):
                raise TypeError("TypeError exception thrown. Expected float in [0,1] for uc_threshold")
            if not (0.0 <= float(uc_threshold) <= 1.0):
                raise ValueError("ValueError exception thrown. Expected uc_threshold to be between 0 and 1")

        if uc_top_k is not None:
            if not isinstance(uc_top_k, int):
                raise TypeError("TypeError exception thrown. Expected integer for uc_top_k")
            if uc_top_k < 1:
                raise ValueError("ValueError exception thrown. Expected uc_top_k to be >= 1")

    # method / transformation
    if uc_method not in ("pearson", "spearman"):
        raise ValueError("ValueError exception thrown. Expected uc_method to be 'pearson' or 'spearman'")

    if uc_transformation not in ("raw", "relative_abundance", "clr"):
        raise ValueError("ValueError exception thrown. Expected uc_transformation to be 'raw', 'relative_abundance', or 'clr'")

    # seed
    if seed is not None and not isinstance(seed, int):
        raise TypeError("TypeError exception thrown. Expected int or None for seed")

    # coracle kwargs
    if coracle_kwargs is None:
        coracle_kwargs = {}
    elif not isinstance(coracle_kwargs, dict):
        raise TypeError("TypeError exception thrown. Expected dict or None for coracle_kwargs")

    # check dimensions
    if x.ndim != 2 or tax.ndim != 2:
        raise ValueError("ValueError exception thrown. Expected x and tax to have two dimensions")

    if x.shape[1] != tax.shape[0]:
        raise ValueError("ValueError exception thrown. Expected x and tax to have the same number of features")

    if x.shape[0] != (y.shape[0] if isinstance(y, pd.DataFrame) else y.shape[0]):
        raise ValueError("ValueError exception thrown. Expected x and y to have the same number of samples")
    ###########################################################################
    
    
    ### check order/hierarchy of taxonomy
    ###########################################################################
    # get tax level names
    tax_levels = list(tax.columns)
    
    # merge to restrict uniqueness counts to ASVs present in aSV
    x_tax = x.transpose().merge(tax, how="left", left_index=True, right_index=True)
    
    # unique counts per level (with prints like before)
    nuniq = {}
    for lvl in tax_levels:
        print(lvl)
        nuniq[lvl] = x_tax[lvl].nunique(dropna=True)
        print(nuniq[lvl])
    
    # detect order left→right
    vals = [nuniq[lvl] for lvl in tax_levels]
    strict_asc  = all(vals[i] < vals[i+1] for i in range(len(vals) - 1))  # highest→lowest
    strict_desc = all(vals[i] > vals[i+1] for i in range(len(vals) - 1))  # lowest→highest (e.g., ASV first)
    
    if strict_desc:
        # if columns go lowest→highest, reverse to make them highest→lowest
        tax = tax.loc[:, tax_levels[::-1]].copy()
        print("taxonomic order has been reversed")
    elif strict_asc:
        # already highest→lowest, do nothing
        pass
    else:
        # ambiguous: keep order unchanged, same message as your original
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! "
              "The order will be unchanged and used as if descending from left to right")
    
    # refresh to reflect any reordering
    tax_levels = list(tax.columns)
    ###########################################################################
        
        
    ### UniCor to filter ASV
    ###########################################################################
    if uc == True: #if UniCor is wanted
        filtered_x = unicorp(x, y, tax, threshold=uc_threshold, top_k=uc_top_k,  method=uc_method, transformation=uc_transformation)
    else: #if not
        x = x.loc[:, (x != 0).any(axis=0)]
        filtered_x = x.transpose().merge(tax, how="left", left_index=True, right_index=True)
        
    ###########################################################################
    
    ### Top-Down Coracle propagation
    ###########################################################################
    level_logs = []
    _level_t0 = time.time()
    _remaining_created = 0

    
    for i in range(len(tax_levels)): #for every taxonomic level
        print(tax_levels[i])
        selected = [] #list of selected groups to propagate
        highest_level = filtered_x.groupby([tax_levels[i]]).sum(numeric_only=True).transpose() #start at the highest level
        highest_level = y.merge(highest_level, left_index=True, right_index=True) #merge with target variable
        x = highest_level.iloc[:,1:] 
        y = highest_level.iloc[:,0].to_frame()
        
        
        resultO = coracle.coracle(x, y, random_state=seed, **coracle_kwargs) #run coracle
        result = resultO.iloc[3:,0] #select only the groups and their respective score
        
        # align scores with the current parent columns
        _scores = result.reindex(x.columns).fillna(0.0)
        
        # bookkeeping: parent counts
        _parents_evaluated = len(x.columns)
        _parents_positive  = int((_scores > 0.0).sum())

        
        if result.max() <= 0:
            print(f"[ZeroScore] halt at '{i}'")
            break
        if i == len(tax_levels)-1:
            print(i)
            break
        # Loop through the top performing groups
        for j in range(len(result)):
            print(result.index[j], result[j])
            if len(selected) >= n_features: 
                # If we already have enough features, break
                break
            if result[j] > 0:
                # If score is not null
                group_features = filtered_x[filtered_x[tax_levels[i]] == result.index[j]]
                group_features = group_features.groupby([tax_levels[i+1]]).sum(numeric_only=True)
                # Get features for this group
                if len(selected) + len(group_features) <= n_features-1:  # If adding all features won't exceed limit
                    selected.extend(group_features.index.tolist())  # Add all features to selected
                    print("adding", group_features.index.tolist())
                else:
                    # Sort remaining features by abundance
                    group_features = group_features.sum(axis=1, numeric_only=True).sort_values(ascending=False)
                    n_to_add = n_features - len(selected)-1 #-1 to leave room for the accumulated features
                    selected.extend(group_features.index[:n_to_add].tolist())
                    print("abundance adding", group_features.index[:n_to_add].tolist())# Add top features by abundance
                    break  # Stop once we reach the required number of features
            else:
                if not selected:
                    raise IndexError("IndexError exception thrown. No score above zero detected...")
                else:
                    print("Warning: threshold partition will not be fulfilled...")
        

        # Accumulate remaining features to ensure compositionality
        if len(selected) == n_features-1:
            filtered_x.loc[~filtered_x[tax_levels[i+1]].isin(selected), tax_levels[i+1]] = "remaining features"
            selected.extend(["remaining features"])
            _remaining_created += 1
            
            
            if len(selected) < n_features:
                print("Warning: final selection has fewer than", n_features, "features: ", len(selected))
        
        print(len(selected), selected)
        filtered_x = filtered_x[filtered_x[tax_levels[i+1]].isin(selected)] #propagate them to the next level
        
        level_logs.append({
            "parent_level": tax_levels[i],
            "child_level": tax_levels[i+1] if i < len(tax_levels) - 1 else None,
            "parents_evaluated": _parents_evaluated,
            "parents_positive": _parents_positive,
            "number_of_children_selected": len(selected),
            "remaining_created": _remaining_created,
            "children_selected": selected,
            "n_features_cap": n_features,
            "elapsed_s": round(time.time() - _level_t0, 4),
        })
        
    ###########################################################################
    unicoracle_log = pd.DataFrame(level_logs)
    print("\nTDS per-level log:\n", unicoracle_log)


    return resultO#return final coracle result


