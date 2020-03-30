import itertools, sys, os
import csv
import glob
import numpy as np
import pandas as pd
import statistics
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
from itertools import tee
import collections

####define scoring scheme and important variables
id_list = ['V','X','Y','Z']
score_list = [10,5,-3,0]

keys = id_list
values = score_list
domain_dictionary = dict(zip(keys,values))

threshold = 0
window = 5000
sliiiiiide_to_the_right = 50

#Let's loop!
file1 = sys.argv[1]
with open(file1, 'r') as file:
    print('Running file: '+ file.name)
    count=0
    count_start=-sliiiiiide_to_the_right

    x_temp = list(file.read())
    x = x_temp[:-1] #delete "/n" that's at the end of the list because we read in the file not explicitly the line
    total_len = int(len(x))

    #resume!
    letter_list = [domain_dictionary[k] for k in x] #convert to scores
    seq_score_nope = sum(letter_list)

    #avg_score = mean(letter_list)
    blocks = int(((len(x) - window) / sliiiiiide_to_the_right) + 1)
    blocks_2 = blocks + 2 #need this otherwise the last (incomplete/little) block will be cut off!
    #print("you will have " + str(blocks_2) + " windows")

    cols = ['Window', 'Position start', 'Position stop','Pass/Fail', 'Score', 'V_count', 'X_count', 'Z_count', 'Y_count']
    dat = pd.DataFrame(columns = cols)
    #
    for i in range(0, blocks_2 * sliiiiiide_to_the_right, sliiiiiide_to_the_right):
        score_result = sum(letter_list[i:i+window])
        new_let_list = x[i:i+window]

        if score_result >= 0 :
            PF_result = "pass"
        else:
            PF_result = "fail"

        #counts for later
        V_count = new_let_list.count('V')
        X_count = new_let_list.count('X')
        Z_count = new_let_list.count('Z')
        Y_count = new_let_list.count('Y')

        #vars for count columns
        count = count +1
        count_start += sliiiiiide_to_the_right #same as c_s = c_s + siiii...
        count_stop = count_start+window

        #dat.index.name = 'Window'
        #let's plot things!
        dat = dat.append({'Window': count,'Position start' : count_start, 'Position stop': count_stop,'Pass/Fail': PF_result, 'Score': score_result, 'V_count': V_count,
                          'X_count': X_count, 'Y_count': Y_count, 'Z_count': Z_count},ignore_index=True)
        #dat.index.name = 'Window'

        outname = (str(file.name)+".tableout.tsv")
        #dat.to_csv(outname, sep='\t', index=False)

        #FIGURES
        pdf_outname = (str(file.name)+".figures.pdf")
        #Character ocunts plot
        #figures = PdfPages(pdf_outname)

    dat.to_csv(outname, sep='\t', index=False)
    #MAIN DATAFRAME CREATED, STORED IN DAT

    #Now let's make the smoothed plot
    df_0 = dat

    #median_0 = df_0['Annotation'].median()
    x = df_0['Window']
    y = np.array(df_0['Score'])
    l = df_0['Window'].count()
    df_empty = pd.DataFrame(index=range(l),columns=range(1))
    for col in df_empty.columns:
        df_empty[col].values[:] = 0

    zero=df_empty[0]

    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    #smooth_val == box_plts
    smooth_val = 100 #####we can change this if we want!


    #statement for handling short sequences (error called if len(y) < smoothing value)
    if len(y) <= smooth_val:
        smooth_val = (0.5 * len(y))
    else:
        smooth_val = smooth_val

    smoth = smooth(y,smooth_val)
    idx = np.argwhere(np.diff(np.sign(zero - smoth))).flatten()
    df = pd.DataFrame(zero[idx])
    df = df.reset_index()
    #we will save to figures, but first we need to do the validation steps

    #This is for validating if region is + or -
    df.loc[-1] = 1  # adding a row for first position
    df.index = df.index + 1  # shifting index
    df = df.sort_index()

    #df.iloc[-1] = len(y)

    #last position as last row
    #print(df['index'])
    df.sort_values(by=['index']) #need to sort first otherwise +1 belwo will break things
    new_list = pd.DataFrame(df['index'] + 1) #df['index'][:-1] + 1 #add +1 to all for next position is +/-, except for last position, will throw erre - so it deletes it, we'll add it in later
    #print(new_list)

    #the_val_to_add = df.iloc[-1] - 1

    #new_list = new_list.append(df.iloc[-1] - 1) #beacuse of +1 transformation few lines above
    new_list_2 = new_list['index']

    #new_list = new_list.append(last_val_to_append, ignore_index=True)

    new_y_val = list(smoth[new_list])   #find position y on smooth line


    #assigning pos / neg for that +1 position
    pos_neg_results = []
    for i in new_y_val:
        if i > 0:
            result = '+'
        else:
            result = '-'
        pos_neg_results.append(result)

    #pos_neg_results.append('N/A') #the last value needs this - not anymore

    #print(pos_neg_results)
    #creating dataframe for next steps
    df.drop(df.columns[len(df.columns)-1], axis=1, inplace=True) #to delete last column, unnamed so tricky to get rid of (?) this does it tho
    df['+/- to the right'] = pos_neg_results
    #print(df['+/- to the right'])
    #append +/- and start stop coords from original table
    df.rename(columns={'index': 'Window'}, inplace=True)

    df['Window']=df['Window'].astype(int)
    df_0['Window']=df_0['Window'].astype(int)

    merged_df = df.merge(df_0, how = 'inner', on = ['Window'])

    merged_df = merged_df.drop(['Pass/Fail','Score','V_count','X_count','Z_count','Y_count'], axis = 1)
    merged_df['Chunk_end'] = 'none'
    merged_df['Window midpoint'] = merged_df.iloc[:,[2,3]].median(axis=1)
    merged_df['Window midpoint'] = merged_df['Window midpoint'].astype(int)

    #df edits to accomodate this:
    #we are duplicating the last row of the df to handle a trailing + chunk (w/ no y=0 intercept to close the chunk)
    merged_df = merged_df.append(merged_df[-1:])
    #now need to make it read actual last stop position (this os not rounded per window like the other coords)
    merged_df = merged_df.replace(merged_df.iloc[-1][3],(total_len+1))
    print(merged_df)

    #now let's get the coordinates for the > 0 'chunks'

    #iterate over for true hit testing
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


    #file name to be used later
    actual_file_name_temp = str(file.name[:-17])


    #this is to define the chunks, accounting for all the ways the graph can look
    #note: leading and trailing here mean a chunk at the start or end of the graph that
    ddf_list = []

    for (i1, row1), (i2, row2) in pairwise(merged_df.iterrows()):
        #for a leading chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] == 0 and \
            row1["Position stop"] != (total_len + 1):
                ddf = ["Chunk_" + str(i1), row1["Position start"], row2["Window midpoint"]]
                ddf_list.append(ddf)
        #for a contained chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] != 0 and \
            row1["Position stop"] != (total_len + 1):
                ddf = ["Chunk_" + str(i1), row1["Window midpoint"], row2["Window midpoint"]]
                ddf_list.append(ddf)
        #3. for a trailing chunk
        if row1['+/- to the right'] == '+' and \
            row1["Position start"] != 0 and \
            row1["Position stop"] == (total_len + 1): #old = merged_df.iloc[0,3]
                ddf = ["Chunk_" + str(i1), row1["Window midpoint"], row2["Position stop"]]
                ddf_list.append(ddf)

        #4. for graphs with no leading and no trailing chunk (for graphs with no y = 0 intercept -> this is is
        #a differently-defined statemnt below b/c the empty file gets appended w/ stuff above from older files when
        #it's in the loop, ALSO the criterion gets fulfilled by contained cunks which means duplicate csv rows for chunks (defined diffrently to specifiy the rules)
        if merged_df.iloc[0,1] == '+' and \
            merged_df.iloc[0,2] == 0 and \
            merged_df.iloc[0,3] == (total_len + 1): #if first column last(2nd row) == last -1 then its one chunk
                rep_list = [('Chunk_0', '0', (total_len+1))]
                ddf_list = rep_list
        else:
                ddf_list = ddf_list

    #print(merged_df)
    #print(ddf_list)

    #make chunk csv
    df = pd.DataFrame(ddf_list)
    this_name = str(file.name+"_chunk_coordinates.csv") #used to be fna_name
    df.to_csv(this_name, index = False)

    ###Find optimal location on plot to place validation marker
    #read in virus table

    #file_name_just_stem = file.name[:-4]

    vir_bait_table = str(actual_file_name_temp+'.VIRUS_BAIT_TABLE.txt')
    with open(vir_bait_table, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        lines = list(reader)

    vir_bait_table = pd.DataFrame(lines)
    vir_bait_table['median'] = round(vir_bait_table[[1,2]].median(axis=1))
    vir_bait_table_med_list = list(vir_bait_table['median'])
    #print(vir_bait_table_med_list)

    points_list = []
    for item in vir_bait_table_med_list:
        eq = round(((item - 2500) + 50) / 50)
        if eq >= len(x):
            plot_point = (len(x) - 1) #1 because it can't = len, has to be less
        else:
            plot_point = eq
        #plot_point = round(((item - 2500) + 50) / 50) #this must stay at = window length (not half like we had talked about, it makes illogical values...basically if the coordinate is towards the end, applying a window 'inbetween' can be out of bounds)
        points_list.append(plot_point)

    new_points_list = [1 if i <=0 else i for i in points_list]

    #print(points_list) #each item represents/is the best/closet window that captures the viral hallmark region

    zero=df_empty[0]

    figures = PdfPages(pdf_outname)

    x2 = (points_list)
    plt.plot(x, y, 'o', ms=0.6)
    plt.axhline(0, 0, l)
    #plt.plot(x, smooth(y,3), 'r-', lw=2)
    #p = smooth(y,100)
    plt.plot(x, smooth(y,100), 'c', lw=2)
    plt.plot(x, smooth(y,100), 'y', markevery = (new_points_list), ms=11.0, marker = '*')
    plt.title("Viral region calls")
    plt.xlabel('Window')
    plt.ylabel('Score')
    plt.rc('axes', titlesize=6.8)     # fontsize of the axes title
    plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=5)    # legend fontsize
    plt.rc('figure', titlesize=8)  # fontsize of the figure title
    plt.grid(True)
    idx = np.argwhere(np.diff(np.sign(zero - smooth(y,100)))).flatten()
    plt.plot(x[idx], zero[idx],  'ro', ms=5.0)
    #plt.plot(x[idx], zero[idx], markevery= (points_list), ms=9.0, marker = 'X', color = 'y')
    #plt.plot()
    df = pd.DataFrame(zero[idx])
    plt.plot()
    plt.savefig(figures, format='pdf')
    plt.close()
    df = df.reset_index()
    #print(df)

    mycol = (["#e7ba52", "#637939", "#7b4173", "#d6616b"])
    dat[['V_count','X_count','Y_count','Z_count']].plot(color = mycol) #same = dat.plot(y=['X_count','N_count','R_count','V_count']) , plt.show()        plt.grid(True)
    plt.grid(True)
    plt.xlabel('Window')
    plt.ylabel('Count')
    plt.title('Character counts')
    plt.rc('axes', titlesize=6.8)     # fontsize of the axes title
    plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=5)    # fontsize of the tick labels
    plt.rc('legend', fontsize=5)    # legend fontsize
    plt.rc('figure', titlesize=8)  # fontsize of the figure title
    plt.savefig(figures, format='pdf')
    plt.close()

    figures.close()
