

# Open the binary file in read mode
with open('boid_data.bin', 'rb') as file:
    # Read the entire contents of the file
    data = file.read()
    # You can process the data here
    print(data[0:3])