rule list_subdirectories:
    output:
        f"{sys.path[0]}/../../results/subdirectories.txt"
    run:
        subdirs = [folder for folder in os.listdir("phastest_results") if os.path.isdir(os.path.join("phastest_results", folder)) and folder.startswith("ZZ_")]
        with open(output[0], "w") as f:
            f.write("\n".join(subdirs))