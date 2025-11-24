# Chord-Progression-Filling-using-YKSP-dijkstra
create any chord progression on any key without professional knowledges - prototype
focusing on chord progression itself on diatonic scale - not blues scale, modes, voicing(altering, extensions)

uses Graph structure consisted of 7 nodes(diatonic scale) and Yen's K-Shortest Paths algorithm to find 10 best chord progressions
user picks key (maj/min), harmony rule type (jazz / classical)
user can input like "C ? G ?" then it'll fill like "C F G C"
(actually, it takes integer values - I is 0, VII is 11, ii is 14, ...vii is 23 and ? is '24')

uses hardcoded weight matrix (maj, min, jazz maj, jazz min)
  i created them considering Cadence rules(VI->V, V->I, vii->I and such), tritone substitution, secondary dominant, ii-V-I
  it'll be better if we make it dynamic - not hardcoded weight matrix
  and ofc just use deep learning algorithm and datas from kaggles to count in diversity

in short, a program for myself when my first and third chord sound cool but second sucks
  
