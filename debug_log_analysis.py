import argparse
import glob
import matplotlib.pyplot as plt



def parse_log_line(raw_line):
    line = raw_line.strip().split(":")
    for i, v in enumerate(line):
        if v.startswith("DEBUG"):
            return line[i], line[i+1], line[i+2:]


def main(args):

    input_logs = glob.glob(args.input)
    if len(input_logs) == 0:
        print("Could not find files matching {}".format(args.input))

    # filesize to mem usage
    input_key = 'DEBUG_INPUT_FILESIZE'
    output_key = 'DEBUG_MAX_MEM'
    is_desired_input = lambda key, extra: key == input_key and len(extra) > 0 and extra[0].endswith('bam')
    is_desired_output = lambda key, extra: key == output_key

    # requested mem to used
    # input_key = 'DEBUG_MAX_MEM'
    # output_key = 'DEBUG_JOB_DISK'
    # is_desired_input = lambda key, extra: key == input_key
    # is_desired_output = lambda key, extra: key == output_key


    for log in input_logs:
        xs = []
        ys = []
        with open(log) as input:
            for line in input:
                key, value, extra = parse_log_line(line)
                if is_desired_input(key, extra):
                    xs.append(int(value))
                elif is_desired_output(key, extra):
                    ys.append(int(value))

        assert len(xs) == len(ys)
        plt.scatter(xs, ys, label=log)

    plt.legend()
    plt.xlabel(input_key)
    plt.ylabel(output_key)
    plt.title("{} to {}".format(input_key, output_key))
    plt.show()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', dest='input', required=True, type=str, help="input file (glob)")

    args = parser.parse_args()
    main(args)