
def run():
    for file in input.fasta_files:
        folder_name = os.path.basename(os.path.dirname(file))
        fasta_number = int(re.search(r"(\d+)\.fasta$", file).group(1))
        recognizer_directory = f"phastest_results/{folder_name}/recognizer_output_{fasta_number}"
        if os.path.isdir(os.path.join(recognizer_directory)):
            print(f"{recognizer_directory} exists for {file}")
        else:
            print("Run recognizer with", file, "from ", folder_name)
            command = f'recognizer -f {file} -o recognizer_output_{fasta_number} -rd /home/fernanda/resources_directory'
            result = Popen(command, shell=True).communicate()
            print(result)