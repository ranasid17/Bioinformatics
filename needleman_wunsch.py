import numpy as np


class needleman_wunsch:

    def __init__(self, X, Y, match, mismatch, gap):
        self.sequence1 = X
        self.sequence2 = Y
        self.match_reward = match
        self.mismatch_penalty = mismatch
        self.gap_penalty = gap

    def build_table(self):
        print("Sequence 1:", self.sequence1)
        print("Sequence 2:", self.sequence2)

        # Pt I: initialize (M+1) x (N+1) table
        opt = np.zeros((len(self.sequence1) + 1, len(self.sequence2) + 1))

        # Pt II: initialize first row/col of alignment table
        for i in range(1, len(opt[0])):  # initialize first row
            opt[0, i] = opt[0, i - 1] - 1  # assign progressive gap penalties

        for i in range(1, len(opt)):  # initialize first col
            opt[i, 0] = opt[i - 1, 0] - 1  # assign progressive gap penalties

        # Pt III: calculate optimal score in each cell
        for i in range(1, len(opt)):  # iter across rows
            for j in range(1, len(opt[0])):  # iter across cols

                vertical_score = opt[i - 1, j] + self.gap_penalty  # calc vert score
                horizontal_score = opt[i, j - 1] + self.gap_penalty  # calc horizontal score

                if self.sequence1[i - 1] == self.sequence2[j - 1]:  # check for match
                    diagonal_score = opt[i - 1, j - 1] + self.match_reward  # diagonal score
                else:  # else is mismatch
                    diagonal_score = opt[i - 1, j - 1] + self.mismatch_penalty  # diagonal score

                # select highest (optimal) score from curr iteration
                optimal_score = max([vertical_score, horizontal_score, diagonal_score])
                # append highest score to [i,j]th mtx elt
                opt[i, j] = optimal_score

        opt_score = opt[-1, -1]  # store optimal alignment score (final cell of mtx)

        return opt, opt_score

    def trace_back(self, table):
        first = ''  # init empty string for alignment of X
        second = ''  # init empty string for alignment of Y

        x = len(self.sequence1)  # store len of first seq
        y = len(self.sequence2)  # store len of second seq

        while (x > 0) or (y > 0):
            # Pt I: calculate vert/horizontal/diagonal scores
            curr_score = table[x, y]  # store curr score from curr iteration

            vertical_score = table[x, y - 1] + self.gap_penalty  # calc vert score
            horizontal_score = table[x - 1, y] + self.gap_penalty  # calc horizontal score

            if self.sequence1[x - 1] == self.sequence2[y - 1]:  # check for match
                diagonal_score = table[x - 1, y - 1] + self.match_reward  # diagonal score
            else:  # else is mismatch
                diagonal_score = table[x - 1, y - 1] + self.mismatch_penalty  # diagonal score

            # Pt II: compare calculated scores to curr_score
            if (y > 0) and (curr_score == vertical_score):  # vert path
                first = '-' + first  # update alignment 1
                second = self.sequence2[y - 1] + second  # update alignment 2
                x = x  # update x idx
                y = y - 1  # update y idx

            elif (x > 0) and (curr_score == horizontal_score):  # horizontal path
                first = self.sequence1[x - 1] + first
                second = '-' + second
                x = x - 1  # update x idx
                y = y  # update y idx

            elif (x > 0) and (y > 0) and (curr_score == diagonal_score):  # diagonal path
                first = self.sequence1[x - 1] + first  # update alignment 1
                second = self.sequence2[y - 1] + second  # update alignment 2
                x = x - 1  # update x idx
                y = y - 1  # update y idx

        opt_alignment = np.array([[first], [second]])

        return opt_alignment
