rows_dict = {}
#generating a dict of rows ord() returns the unicode value of a character, chr() returns the character of a unicode value
for i in range(ord('A'), ord('I')):
    row = chr(i)
    rows_dict[row] = []
    
new_dict = {}
for _ in range(ord('A'), ord('I')):
    row = chr(_)
    new_dict[row] = []
    
print(new_dict)