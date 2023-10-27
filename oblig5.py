b = 1/3
k = 1/10

s = 1
i = 10/(5e6)
r = 0

def updateSIR(s_y, i_y, r_y, b, k):
    s_new = s_y - b*s_y*i_y
    i_new = i_y + b*s_y*i_y - k*i_y
    r_new = r_y + k*i_y
    return s_new, i_new, r_new

with open("data/pan2.csv","w") as output:
    output.write("day,susceptible,infected,recovered\n")
    s_in, i_in, r_in = s, i, r
    for i in range(120):
        s_in, i_in, r_in = updateSIR(s_in, i_in, r_in, b, k)
        print(s_in, i_in, r_in)
        output.write(f"{i+1},{s_in:.4e},{i_in:.4e},{r_in:.4e}\n")
