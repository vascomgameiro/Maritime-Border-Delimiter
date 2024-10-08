from src.fileHandler import FileHandler


def main():
    file_handler = None

    while True:
        command = input("?> ").strip().split()
        if len(command) == 0:
            continue
        cmd = command[0]

        if cmd == "import_points":
            if len(command) < 2:
                print("Error: missing file name.")
                continue
            file_name = command[1]

            try:
                with open(file_name, "r"):
                    pass
            except FileNotFoundError:
                print(f"Error: File '{file_name}' does not exist.")
                continue
            except IOError:
                print(f"Error: File '{file_name}' is not accessible or readable.")
                continue

            try:
                file_handler = FileHandler(file_name)
                print(f"Points from {file_name} imported successfully.")
            except Exception as e:
                print(f"Error reading points: {e}")

        elif cmd == "make_map":
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue
            try:
                file_handler.make_map()
                print("Map created successfully.")
            except Exception as e:
                print(f"Error creating map: {e}")

        elif cmd == "quit":
            print("Exiting the program.")
            break

        else:
            print(f"Unknown command: {cmd}")


if __name__ == "__main__":
    main()
