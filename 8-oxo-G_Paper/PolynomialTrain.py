import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

def midpoint_tf_analysis_polyfit(input_file, tf_map):
    sns.set_style("whitegrid")  # Set seaborn style for the plot

    input_path, tf_path = Path(input_file), Path(tf_map)
    output = tf_path.with_name(input_path.stem + '_' + tf_path.stem + 'intersect.txt')
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {output}'],
                        stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    with open(output) as f:
        position_dict = {}
        for line in f:
            tsv = line.strip().split('\t')
            mut_start = int(tsv[1])
            tf_start = int(tsv[9])
            tf_end = int(tsv[10])
            mid_pt = math.floor((tf_end + tf_start) / 2)
            mut_position = mut_start - mid_pt
            position_dict[mut_position] = position_dict.setdefault(mut_position, 0) + 1

        # Sort the data by position
        sorted_data = sorted(position_dict.items())
        positions, counts = zip(*sorted_data)
        positions = np.array(positions)
        counts = np.array(counts)

        # Perform k-fold cross-validation for different polynomial degrees
        k = 10  # Number of folds in cross-validation
        max_degree = 50  # Maximum degree to consider
        degrees = list(range(1, max_degree + 1))
        cv_scores = []

        for degree in degrees:
            model = make_pipeline(PolynomialFeatures(degree), LinearRegression())
            scores = cross_val_score(model, positions.reshape(-1, 1), counts, cv=k)
            cv_scores.append(np.mean(scores))

        # Find the optimal degree
        optimal_degree = degrees[np.argmax(cv_scores)]

        # Train the model with the optimal degree using all the data
        optimal_model = make_pipeline(PolynomialFeatures(optimal_degree), LinearRegression())
        optimal_model.fit(positions.reshape(-1, 1), counts)

        # Plot the data and the fitted curve
        plt.scatter(positions, counts, label="Original data", alpha=0.5)
        plt.plot(positions, optimal_model.predict(positions.reshape(-1, 1)), color="red", linewidth=3, label="Optimal fit")
        plt.xlabel("Position")
        plt.ylabel("Counts")
        plt.legend()
        plt.title(f"Comparison of Original Data and Optimal Polynomial Fit (Degree {optimal_degree})")
        plt.show()

# Call the function with your input_file and tf_map
mutations = '/media/cam/Data9/CortezAnalysis/Cam_calls/Analysis/vcf_files/concat/KM_treated.bed'
old_map = '/media/cam/Data9/CortezAnalysis/Cam_calls/Analysis/TF_Binding/new_tf_map/non_filtered_TF_map_1000.bed'
midpoint_tf_analysis_polyfit(mutations, old_map)
