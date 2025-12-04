"""Execution file."""

import argparse
import src.training as training
import src.plotting as plotting
import src.scoring as scoring

def main():
    parser = argparse.ArgumentParser(description="Run training, plotting, or scoring steps.")

    # Boolean flags for the three tasks
    parser.add_argument("--train", action="store_true", default=False,
                        help="Run training (default: False)")

    parser.add_argument("--plot", dest="plot", action="store_true", default=True,
                        help="Run plotting (default: True)")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="Disable plotting")

    parser.add_argument("--score", dest="score", action="store_true", default=True,
                        help="Run scoring (default: True)")
    parser.add_argument("--no-score", dest="score", action="store_false",
                        help="Disable scoring")

    args = parser.parse_args()

    # Execute based on flags
    if args.train:
        training.run_train()

    if args.plot:
        plotting.make_plot()

    if args.score:
        scoring.run_score()




if __name__ == "__main__":
    main()

