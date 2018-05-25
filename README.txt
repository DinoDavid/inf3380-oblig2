Dette er hva jeg har gjort:

Jeg byttet om på 0 og 1 i Cart_shift

Jeg slettet nullstilling av c i matrix_mult

Jeg fikset på deallocate_matrix funksjonen

Jeg fant ut at jeg manglet large_matrix_b.bin så jeg kopierte den i mappa

Flyttet på moff og noff

Lagde en C_part_m_max og C_part_n_max

Når jeg kalkulerer matrix_max så har jeg med +rest istedenfor +1 som jeg hadde tidligere

Byttet på rekkefølgene på MPI_Sendrecv_replace og fikset litt på dem

Fikset duplikat av col_displ_a

Lagde openMP delen i matrix_mult

Testet om openMP delen fungerer ved å printe ut med en for loop


Koden har vært testet med et positivt resultat for følgende antall processeringsenheter:

small_matrix - 4, 9, 25

large_matrix - 4, 9, 16, 25, 36

Jeg brukte programmet til kompisen min compare.c til å sjekke om matrisene er like og fikk "comparison passed" med alle de over.
