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

s, i ,r = updateSIR(s, i, r, b, k)

data = {"day":[0], "s":[s], "i":[i], "r":[r]}

for i in range(1,121):
    s_new, i_new, r_new = updateSIR(data["s"][-1], data["i"][-1], data["r"][-1], b, k)
    data["day"].append(i+1)
    data["s"].append(s_new)
    data["i"].append(i_new)
    data["r"].append(r_new)

with open("data/pan2.csv","w") as outfile:
    outfile.write("day,susceptible,infected,recovered\n")
    for i in range(120):
        outfile.write(f"{data['day'][i]},{data['s'][i]:.4e},{data['i'][i]:.4e},{data['r'][i]:.4e}\n")
