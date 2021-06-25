import matplotlib.pyplot as plt

pie_data = {'0' : 37,
           '1' : 5,
           '2' : 2,
           '?' : 10}

bar_data = {23 : 3,
            29 : 2,
           31 : 2,
           41 : 1,
           47 : 1,
           73 : 1
           }

def plot_pie():
    labels = []
    sizes = []

    for x, y in pie_data.items():
        labels.append(x)
        sizes.append(y)

    plt.pie(sizes, labels=labels)

    plt.axis('equal')
    plt.title('Number of new isogeny primes')
    plt.show()


def plot_bar():

    plt.bar(*zip(*bar_data.items()))
    # plt.axis('equal')
    plt.xticks(list(bar_data.keys()))
    plt.title('Frequency of new isogeny primes')
    plt.show()

if __name__ == "__main__":
    plot_pie()
    plot_bar()