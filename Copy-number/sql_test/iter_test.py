query = (('K562','NAT2',3),('K562','NMNAT2',4),('K562','GNAT2',4))

hit_num = len(query)
inside_index=0
i = 0
#number of elements in 'hit' tuple, init @ 0
print(f"Number of hits: {hit_num}")
hit=[]

    
while i < hit_num:
    while inside_index < 3:
        hit.append(query[i][inside_index])
        inside_index += 1
    print(hit)
    hit=[]
    inside_index = 0
    i += 1
        

    

