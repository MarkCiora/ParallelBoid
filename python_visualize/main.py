import struct

# Open the binary file in read mode
with open('boid_data.bin', 'rb') as file:
    # Read the entire contents of the file
    data = file.read()
    # You can process the data here
    print(data[1:5])
    print(struct.unpack('i', data[0:4]))
    print(struct.unpack('i', data[4:8]))
    print(struct.unpack('f', data[8:12]))
    print(struct.unpack('f', data[12:16]))