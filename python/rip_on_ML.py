from ripser import ripser

def main(d_matrix):
    bars = ripser(d_matrix,  distance_matrix=True, maxdim=1)
    PD = bars['dgms']
    return(PD)