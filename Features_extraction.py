data_p1 = pd.read_csv("aa03l_preprocessed_eyedata_SWITCH.csv", sep = ";")
data_p2 = pd.read_csv("AA30W_preprocessed_eyedata_SWITCH.csv", sep = ";")

def feature_generation(data):
    Timestp = data.groupby("TrialNr")["Time"].min()

    Total_trial_blinks = data.groupby("TrialNr")["eye_event"].apply(lambda x :
                            (x == "start_blink").count())
    #number of blinks

    Average_pup_size = data.groupby("TrialNr")["pupil"].mean()

    SD_pupil = data.groupby("TrialNr")["pupil"].std()

    SD_gaze_x = data.groupby("TrialNr")["gazeX"].std()
    SD_gaze_y = data.groupby("TrialNr")["gazeY"].std()
    #to see how much the gaze deviates
    #we could maybe do gaze average as well?

    Fix_disp_proportion = data.groupby("TrialNr")["fixdisp"].apply(lambda x : ((x < 0).sum()) /x.count() )
    #propotion of crossed eyes
    Blink_proportion = data.groupby("TrialNr")["eye_cat"].apply(lambda x : ((x =="blink").sum()) /x.count() )
    #proportion of the recordings that were within a blink
    Fixation_proportion = data.groupby("TrialNr")["eye_cat"].apply(lambda x : ((x =="fixation").sum()) /x.count() )
    #proportion of the recordings that were within a fixation
    Saccade_proportion = data.groupby("TrialNr")["eye_cat"].apply(lambda x : ((x =="saccade").sum()) /x.count() )
    #proportion of the recordings that were within a saccade

     #not sure if we should count the NAs while doing proportion,
    #if we do replace x.count() by len(x)

    Microsacc_proportion = data.groupby("TrialNr")["micros"].apply(lambda x : ((x ==1).sum()) /len(x) )
    #proportion of the recordings that were within a microsaccade

    #Trialpos = data.groupby([('block'=="block1"), 'TrialNr'])["TrialPos"].unique()
    #Trial_block = data["trial_block"]
    #Trial = data["TrialNr"]
    #Trials are not ordered (see timestamp) and 2 same blocks in all trials

    features = np.column_stack([
        Timestp, Average_pup_size, SD_pupil, SD_gaze_x, SD_gaze_y, Total_trial_blinks,Blink_proportion, Fix_disp_proportion,  Fixation_proportion, Saccade_proportion, Microsacc_proportion
    ])

    features_names = ["Time_trial_begins", "Average_pupil_size", "SD_pupil", "SD_gaze_x", "SD_gaze_y", "Total_trial_blinks","Blinking_proportion", "Cross_eyes_proportion",  "Fixation_proportion", "Saccade_proportion", "Microsacc_proportion"]
    return features, features_names

features_p1, feature_names = feature_generation(data_p1)
features_p2, _ = feature_generation(data_p2)


