import os

files = os.listdir()

for f in files:
    f_split = f.split('_')

    try:
        if int(f[0]):
            f_end = "-1.mat"
            f_split[-1] = f_end
            f_unsplit = ""
            for i, fs in enumerate(f_split):
                f_unsplit += fs + ("_" if i < (len(f_split)-1) else "")
            print(f_unsplit)
            os.rename(f, f_unsplit)
    except ValueError:
        print(f_split[0], " is not a number")
