
def hamming_distance(s1, s2):
    a = 0
    
    if len(s1) == len(s2):
        for index in range (0,len(s1)):
            if s1[index] != s2[index]:
                a+=1
        return a
                
    else:
        print("ValueError")


if __name__ == "__main__":
    s1 = input("염기서열을 입력하시오.")
    s2 = input("염기서열을 입력하시오.")
    print(hamming_distance(s1, s2))
