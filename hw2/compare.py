# Exercise 6.1. Write a compare function that returns 1 if x > y, 0 if x == y, and -1 if x < y.

def compare(x, y):
    if x > y:
        return 1
    elif x == y:
        return 0
    elif x < y:
        return -1

if __name__ == "__main__":
    print compare(4, 8)