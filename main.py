import argparse
from math import ceil
import src.training as training
import src.plotting as plotting
import src.scoring as scoring
import utils.model as model


def main():
    parser = argparse.ArgumentParser(description="RNA scoring pipeline")

    # ===== FLAGS: all steps run by default =====
    parser.add_argument("--no-train", action="store_true",
                        help="Disable the training step")
    parser.add_argument("--no-plot", action="store_true",
                        help="Disable the plotting step")
    parser.add_argument("--no-score", action="store_true",
                        help="Disable the scoring step")

    # ===== DIRECTORIES =====
    parser.add_argument("--trainset", default="data/structures/train",
                        help="Directory containing training PDB/CIF files")
    parser.add_argument("--profiles", default="data/profiles",
                        help="Directory to save or load trained profiles")
    parser.add_argument("--testset", default="data/structures/test",
                        help="Directory containing test PDB/CIF files")
    parser.add_argument("--scores", default="data/scores",
                        help="Directory to store scoring results")

    # ===== MODEL PARAMETER OVERRIDES =====
    parser.add_argument("--max-distance", type=int, default=None,
                        help="Set maximum allowed distance cutoff (default: 20)")
    parser.add_argument("--position-skip", type=int, default=None,
                        help="Minimum residue separation (default: 4)")
    parser.add_argument("--maximum-score", type=int, default=None,
                        help="Maximum allowed score value in the statistical potential (default: 10)")
    parser.add_argument("--bin-width", type=float, default=None,
                        help="Histogram bin width for distance distributions (default: 1.0 Å)")

    args = parser.parse_args()

    # ===== APPLY MODEL PARAMETER OVERRIDES =====
    if args.max_distance is not None:
        model.max_distance = args.max_distance
        model.max_distance_sq = args.max_distance ** 2

    if args.position_skip is not None:
        model.position_skip = args.position_skip

    if args.maximum_score is not None:
        model.maximum_score = args.maximum_score

    if args.bin_width is not None:
        model.bin_width = args.bin_width
        model.num_bins = ceil(model.max_distance / model.bin_width)  

    # ===== DETERMINE WHICH STEPS SHOULD RUN =====
    run_training = not args.no_train
    run_plotting = not args.no_plot
    run_scoring = not args.no_score

    # ===== PRINT PIPELINE SUMMARY =====
    print("\n=== PIPELINE CONFIGURATION ===")
    print("Training step: ", run_training)
    print("Plotting step: ", run_plotting)
    print("Scoring step:  ", run_scoring)
    print("Training set directory:  ", args.trainset)
    print("Profile directory:       ", args.profiles)
    print("Test set directory:      ", args.testset)
    print("Scores output directory: ", args.scores)
    print("================================\n")

    # ===== RUN THE SELECTED STEPS =====
    if run_training:
        training.run_train(args.trainset, args.profiles)

    if run_plotting:
        plotting.make_plot()

    if run_scoring:
        scoring.run_score(args.profiles, args.testset, args.scores)


if __name__ == "__main__":
    main()
