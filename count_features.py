from sys import argv

if __name__ == "__main__":

    with open(argv[1], "r") as ifi:
        for line in ifi:
            if not "Mut" in line:
                tokens = line.strip().split("\t")
                key = tokens[0]
                vals = sum([int(i) for i in tokens[1:]])
                print(key, vals)
