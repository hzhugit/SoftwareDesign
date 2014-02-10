def check_fermat(a, b, c, n):
    if n > 2 and a**n + b**n == c**n:
        print "Holy smokes, Fermat was wrong!"
    else:
        print "No, that doesn't work."

def fermat_input():
    a = int(raw_input("Input a: "))
    b = int(raw_input("Input b: "))
    c = int(raw_input("Input c: "))
    n = int(raw_input("Input n: "))
    check_fermat(a, b, c, n)

fermat_input()