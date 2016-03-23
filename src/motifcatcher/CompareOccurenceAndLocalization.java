package motifcatcher;

/**
 * Created by IntelliJ IDEA.
 * User: cobalt
 * Date: 06.12.2011
 * Time: 13:58
 * To change this template use File | Settings | File Templates.
 */
public class CompareOccurenceAndLocalization {

    public class CoRes {
        String CoOccurrences;
        String CoLocalizations;
        String CoVals;
    }

    public CoRes MC_CompareOccurrenceAndLocalization(MotifMap mm) {

        //This function determines the significant co-occurrence / co-localizations
        //of instances of families of motifs on different sites in the superset.
        //
        //Inputs:
        //   MotifMap
        //       This structure contains all sites with any instances of familial
        //       motifs, as well as the PWMs of those families.
        //   3
        //       Number of initial rows reserved in MotifMap structure for Header.

        //
        //Outputs:
        //   CoOccurrences:
        //       Structure emphasizing Co-Occurrence.  Cases are sorted by their
        //       respective Co-Occurrence scores.
        //   CoLocalizations:
        //       Structure emphasizing Co-Localization.  Cases are sorted by their
        //       respective Co-localization scores.
        //   CoVals:
        //       Emphasizes a combination of Co-Occurrence and Co-Localization.
        //       Cases are sorted by their combined score.
        //
        //Algorithm:
        //
        // 'Significance' here is a combination of [1] co-occurrence and [2]
        // co-localization.  However, considering only one or the other may be more
        // important, so output is also provided entirely on one metric or the
        // other.
        //
        // [1] - Co-occurrence
        //
        // The probability of a single event occurring (instance of a familial motif
        // occurring on a site) is computed by the observed frequency.  Combined
        // events are also counted, generating a long list of all observed cases,
        // and the frequencies.
        //
        // From these figures, probabilities are progressively computed (by solving
        // more and more complicated algebraic equations) determining the
        // probability of a particular co-dep}ence (the combined effect of
        // multiple families occurring).  Note that if familial occurrence is truly
        // indep}ent, combined occurrence probabilities should all be zero.
        //
        // The relative contribution of the probability of multiple families
        // occurring over all possible ways to generate the same result is computed,
        // and ranked.  Families that show a high degree of dep}ence with each
        // other will demonstrate a higher fraction of the combined probabilities of
        // each versus all probabilities.
        //
        //
        // CAUTION:
        //   The deeper conditional probabilities are computed, the more
        //   opportunities there are for error, and for error propagation. Take
        //   highly-nested co-occurrence probabilities with care.
        //
        // CAUTION:
        //   Computing highly-nested co-occurrence probabilities is very
        //   computationally intensive - not recomm}ed to compute beyond layer 6.
        //   this catch is built into the code!!
        //
        // [2] - co-localization
        //
        // The compariatively easier co-localization problem is determined by
        // examining each intersecting subset that has a frequency greater than one,
        // and determining the distance from the center of motifs of each type among
        // all co-occurring families.
        //
        // Pearson's correlation coefficient is determined for each individual
        // co-variance among motifs, and averaged over all co-occurring
        // relationships.  All co-localizations are considered equally weighted.
        //
        // Major Steps:
        //   (0) Testing/Initializations
        //
        //   [1] - Co-Occurrences
        //
        //   (1) Determination of initial (single familial instance) probabilities.
        //   (2) Determiniation of familial instance profiles for all sites in the
        //       superset.
        //   (3) Observed frequency values are added to the list.  Note that this
        //       is the step with greatest uncertainty, and can result in the most
        //       error.
        //   (4) Determine the dep}ent probabilities for all cases that involve
        //       instances of more than one motif family.
        //   (5) Evaluate Co-Dep}ence by the fraction of specific co-dep}ence
        //       versus all possible probabilities to generate a particular state.
        //
        //   [2] - Co-Localizations
        //
        //   (6) Determine Co-Localizations
        //   (7) Evaluate / score / rank Co-Localizations
        //
        //   [3] - Co-Vals
        //
        //   (8) Calculate combined Co-Occurrence and Co-Localization values
        //       (Co-Vals).
        //       CoVals are computed by:
        //
        //       Frequency * (PCC+b1)^2 * (CoDep+b2)
        //           note b1 = 0.2
        //                b2 = 0.5
        //       (subject to change)

        ////
        // ----------------------------------------------------------------------- //
        // (1) Determine probabilities of finding each family.
        // ----------------------------------------------------------------------- //
/*
        //initialize
        FamilyProbs = zeros(1,length(MotifMap(1,:))-1);

        //determine whole set length
        WholeSetLength = length(MotifMap(:,1)) - 3;

        //Counter to keep track of inputs
        int Counter = 0;
        for i = 2:length(MotifMap(1,:))
           NumDetermined = 0;

           for j = 3+1:length(MotifMap(:,1))
                if ~isempty(MotifMap{j,i})
                    NumDetermined = NumDetermined + 1; //increment
                }
           }

           //write results to FamilyProbs structure.

           //probability of discovering a given motif family in a site in the
           //superset
           Counter = Counter + 1;
           FamilyProbs(Counter) = NumDetermined/WholeSetLength;

        }
        ////
        // ----------------------------------------------------------------------- //
        // (2) Determine motifs involved in each site, and number of occurrences.
        // ----------------------------------------------------------------------- //

        MotifsInSites = {'Motifs','Frequency','Observed Frequency',...
            'Expected Frequency (exact co-dep}ence)','Co-Occurrence Fraction'};
        Counter = 1;
        for i = 3+1:length(MotifMap(:,1))
            Motifs = [];
            for j = 2:length(MotifMap(1,:))
                if ~isempty(MotifMap{i,j})
                    //add occurrence to the list of motifs for this set.
                    Motifs = [Motifs j-1];
                }
            }

            if ~isempty(Motifs)

                //segment 'Motifs' into all possible subsets.
                Motifs = Subsets(Motifs);

                //account for each possible subset, and add to output structure.
                for k = 1:length(Motifs(:,1))

                    DoNotAdd = true;
                    if Counter ~= 1
                        for j = 1:length(MotifsInSites(:,1))
                           if isequal(Motifs{k,1},MotifsInSites{j,1})
                               DoNotAdd = false;
                               MotifsInSites{j,2} = MotifsInSites{j,2} + 1;
                           }
                        }
                    }

                    if DoNotAdd == true
                       Counter = Counter + 1;
                       MotifsInSites{Counter,1} = Motifs{k,1};
                       MotifsInSites{Counter,2} = 1;
                    }

                }

            else //case of no motifs observed - so, no subsets.

                    DoNotAdd = true;
                    if Counter ~= 1
                        for j = 1:length(MotifsInSites(:,1))
                           if isequal(Motifs,MotifsInSites{j,1})
                               DoNotAdd = false;
                               MotifsInSites{j,2} = MotifsInSites{j,2} + 1;
                           }
                        }
                    }

                    if DoNotAdd == true
                       Counter = Counter + 1;
                       MotifsInSites{Counter,1} = Motifs;
                       MotifsInSites{Counter,2} = 1;
                    }
            }
        }

        //sort list (bubble sort)
        //first, by length
        for i = 2:Counter-1
           for j = 2:Counter-1
              if length(MotifsInSites{j,1}) < length(MotifsInSites{j+1,1})
                  temp = MotifsInSites(j,:);
                  MotifsInSites(j,:) = MotifsInSites(j+1,:);
                  MotifsInSites(j+1,:) = temp;
              }
           }
        }

        //then, within a length class, by ordering of components.

        ////
        // ----------------------------------------------------------------------- //
        // (3) Add Observed Frequency values
        // ----------------------------------------------------------------------- //
        for i = 2:Counter

            if ~isempty(MotifsInSites{i,1})

               //Compute observed probability
               Freq = MotifsInSites{i,2}/WholeSetLength;

               MotifsInSites{i,3} = Freq;

               //for single cases, the expected probability is the same as the
               //observed frequency (by definition).
               if length(MotifsInSites{i,1})==1
                  MotifsInSites{i,4} = MotifsInSites{i,3};
               }

            else

                //Compute expected probability
                Probability = 1;
                for j = 1:length(FamilyProbs)
                   Probability = Probability * (1-FamilyProbs(j));
                }

                MotifsInSites{i,4} = Probability;

                //Compute observed Probability
                Freq = MotifsInSites{i,2}/WholeSetLength;

                MotifsInSites{i,3} = Freq;
            }
        }
        ////
        // ----------------------------------------------------------------------- //
        // (4) Compute expected probability values for nonempty cases
        // ----------------------------------------------------------------------- //

        //iterate through matrix, working through every N-ple class
        //repeated process for all n-ples > 1

        //comment: computing dep}encies beyond n = 6 takes a very long time
        for j = length(MotifsInSites(:,1)):-1:2 //scan whole set

            if length(MotifsInSites{j,1}) > 1  && length(MotifsInSites{j,1}) < 7
                //evaluate this case.
                //we only evaluate cases of the same length at the same time.
                //the set is ordered, so iterating through accomplishes this.

                 //First, compute all 'True' variables from the list.
                 //initialize
                 TrueVal = Subsets(MotifsInSites{j,1});
                 TrueVal{1,2} = 'x'; //variable to solve for
                 for k = 2:length(TrueVal(:,1))
                 //for k = 2:2+nchoosek(length(MotifsInSites{j,1}),length(MotifsInSites{j,1})-1)-1;

                        //first, get exact match.  Write this term.
                        for z = 2:Counter
                           if isequal(TrueVal{k,1},MotifsInSites{z,1})
                               TrueVal{k,2} = strcat('(',num2str(MotifsInSites{z,4},14));
                           }
                        }

                        //now, retrieve all subsets, from

                        for z = 1:length(TrueVal(:,1))
                           if isequal(intersect(TrueVal{z,1},...
                                   TrueVal{k,1}),TrueVal{k,1}) && ...
                                   ~isequal(TrueVal{z,1},...
                                   TrueVal{k,1})

                               TrueVal{k,2} = strcat(TrueVal{k,2},'-',TrueVal{z,2});
                           }
                        }

                        //close the set.
                        TrueVal{k,2} = strcat(TrueVal{k,2},')');
                 }

                 //determine all possible products that would amount to the
                 //probability of selection, and write this into the
                 //equation.
                 AllProducts = SetPartition(MotifsInSites{j,1});

                 eqn = '';
                 for k = 1:length(AllProducts) //each possible product

                     //initialize eqn term
                     eqnterm = '';

                     //build eqn term
                     for q = 1:length(AllProducts{k,1})

                            for s = 1:length(TrueVal(:,1))
                                if isequal(TrueVal{s,1},AllProducts{k,1}{1,q})
                                   term =  TrueVal{s,2};
                                }
                            }

                            if strmatch(eqnterm,'')
                                eqnterm = term;
                            else
                                eqnterm = strcat(eqnterm,'*',term);
                            }

                     }

                     //enclose each eqnterm in parantheses.
                     eqnterm = strcat('(',eqnterm,')');

                     //add eqnterm to eqn.
                     if ~isempty(eqn)
                        eqn = strcat(eqn,'+',eqnterm);
                     else
                        eqn = eqnterm;
                     }

                 }

                 eqn = strcat(eqn,'=',num2str(MotifsInSites{j,3},14));

                 //once the equation has been assembled, solve it
                 A = solve(eqn);

                 //candidate solutions must be real, greater than 0, and
                 //less than 1. Candidate solutions must also not be larger than
                 //the minimal individual case probability they contain.

                 //Finally, candidate solutions may not be larger than the
                 //observation frequency for that case; as the candidate
                 //solution represents the part of all ways to generate a
                 //particular state that is entirely the simultaneous
                 //co-occurrence of all included states.

                 minval = 1;
                 for k = 1:length(MotifsInSites{j,1})
                     if FamilyProbs(MotifsInSites{j,1}(k)) < minval
                         minval = FamilyProbs(MotifsInSites{j,1}(k));
                     }
                 }


                 sol = [];
                 for k = 1:length(A)
                     if isreal(A(k))
                         if double(A(k))>0 && double(A(k))<minval && double(A(k))<MotifsInSites{j,3}
                            sol = [sol double(A(k))];
                         }
                     }
                 }

                 //write to output structure.
                 if length(sol) > 1

                     //acceptable case:
                     //smallest probability satisfies additional
                     //maximization/minimization constraints of nested
                     //probabilities.
                     sol = min(sol);

                     MotifsInSites{j,4} = sol;
                 elseif length(sol) == 1

                     //ideal case:
                     //only one relevant solution found
                     MotifsInSites{j,4} = sol;
                 else

                     //problem case:
                     //stochastic error in 'observed
                     //frequency' term yields no positive roots; also,
                     //accumulation of errors of this kind from lower-order
                     //solved equations with stochastic errors in 'observed
                     //frequency' terms.
                     //
                     //In that case, we just set the probability equal to
                     //zero (no co-dep}ency).
                     MotifsInSites{j,4} = 0;
                 }

            else

                    //for very high-order nesting, just set the emerging dep}ency
                    //to zero.
                    MotifsInSites{j,4} = 0;

            }

            //optional print statement
            disp(['Class ' num2str(length(MotifsInSites{j,1})) ': Case ' num2str(j) ' Evaluated.'])
        }
        ////
        // ----------------------------------------------------------------------- //
        // (5) Evaluate and sort significant Co-Occurrences
        // ----------------------------------------------------------------------- //


        //Finally, to make sense of the data, compute the degree of occurrence of a
        //particular state based solely on the co-occurrence of all component
        //states, as compared to all observed degrees of occurrence.
        //
        //This value describes the fraction of occurrence where co-occurrence is the
        //key factor.
        //
        //Note that this value is perfect for all 2-ples, usually well-determined
        //for 3-ples, and proceeds with increasing error up to n-ples.  The
        //increased error is due to uncertainty in the observed frequencies for
        //component cases - observed frequencies might not be consistent with one
        //another.
        //
        //Cases where this occurs are not considered to be significant as far as
        //co-occurrence evaluations are concerned.
        for j = 2:Counter
            MotifsInSites{j,5} = (MotifsInSites{j,4}/MotifsInSites{j,3});
        }

        //sort by the most significiant co-occurrences.  Significance of
        //co-occurrence is represented here by the relative importance of the
        //particular evaluated intersection.
        Header = MotifsInSites(1,:);
        Body = sortrows(MotifsInSites(2:Counter,:),5);

        CoOccurrences = vertcat(Header,flipud(Body));
        ////
        // ----------------------------------------------------------------------- //
        // (6) Determine Co-localization among co-occurring cases
        // ----------------------------------------------------------------------- //

        //build a special cell array that examines all co-occurrences for all cases
        //where things co-occur more than once.

        CoLocalizations(1,:) = MotifsInSites(1,:);
        CoLocalizations{1,6} = 'Motif Centers';

        RowCounter = 1;

        for j = 2:Counter
           if CoOccurrences{j,2} > 1 && ~isempty(CoOccurrences{j,1}) && length(CoOccurrences{j,1}) > 1
                RowCounter = RowCounter + 1;
               CoLocalizations(RowCounter,1:5) = CoOccurrences(j,:);
           }
        }

        //for each significant co-occurrence, compute the distance.
        //find the cases in the MotifMap structure.

        for i = 2:RowCounter //repeat procedure for every row.

            //initializations: every row has a case
            Case = CoLocalizations{i,1};

            EntryCounter = 1;

            //Title
            for j = 1:length(Case)
                CoLocalizations{i,6}{1,j} = strcat('Family ',num2str(Case(j)),':');
            }

            for j = 3+1:length(MotifMap(:,1))

                Motifs = []; //initialize for each row.

                //recover motifs from each row.
                for k = 2:length(MotifMap(1,:))
                    if ~isempty(MotifMap{j,k})
                        //add occurrence to the list of motifs for this set.
                        Motifs = [Motifs k-1];
                    }
                }

                //test row to see if that row corresponds to a superset of 'Case'.
                if isequal(intersect(Motifs,Case),Case)

                    //optional print statement
                    //disp(num2str(j))

                    //increment entry counter.
                    EntryCounter = EntryCounter + 1;
                    for k = 1:length(Case)
                        CoLocalizations{i,6}{EntryCounter,k} = ...
                            mean(MotifMap{j,Case(k)+1});
                    }
                }

            }

        }

        ////
        // ----------------------------------------------------------------------- //
        // (7) Evaluate Co-localization among motif centers
        // ----------------------------------------------------------------------- //

        //each CoLocalizations{i,6}(2:length(CoLocalizations{i,6}(:,1)),:)

        CoLocalizations{1,7} = 'Distance Matrix';
        CoLocalizations{1,8} = 'Co-Localization (PCC)';

        for i = 2:RowCounter
            //extract each matrix containing centers
            Centers = (cell2mat(CoLocalizations{i,6}(2:length(CoLocalizations{i,6}(:,1)),:)));

            //initialize the Covariance distance matrix
            CoDistanceMatrix = zeros(length(Centers(:,1)),nchoosek(length(Centers(1,:)),2));

            ColCounter = 0;
            for j = 2:length(Centers(1,:))
               for k = 1:j-1
                   ColCounter = ColCounter + 1;
                   CoDistanceMatrix(:,ColCounter) = (Centers(:,j)-Centers(:,k));
               }
            }

            //critical step: switch rows and columns to put values into correct form
            //for MATLAB's built-in correlation computation function.
            CoDistanceMatrix = CoDistanceMatrix';

            //write to output array.
            CoLocalizations{i,7} = CoDistanceMatrix;


            //calculate Pearson's correlation coefficient for each pair of columns.
            //
            //assume an equal weight amongst all correlations, and determine the
            //average correlation among the whole set.
            //this is achieved by taking the lowerleft triangle of the pxp
            //correlation matrix.

            CoLocalizations{i,8} = mean(abs(LowerLeft(corr(CoDistanceMatrix))));

            //if the Pearson Correlation Coefficient cannot be computed (because of
            //0-row or 0-column, the result is 'NaN', so we consider this case to be
            //unimportant for our analysis.
            if isnan(CoLocalizations{i,8})
               CoLocalizations{i,8} = 0;
            }
        }


        Header = CoLocalizations(1,:);
        Body = sortrows(CoLocalizations(2:RowCounter,:),8);

        CoLocalizations = vertcat(Header,flipud(Body));

        ////
        // ----------------------------------------------------------------------- //
        // (8) Combined CoVal
        // ----------------------------------------------------------------------- //

        //Scores are based on
        //   (1) Frequency number (number of co-occurrences)
        //   (2) Co-Occurrence (computed dep}ency)
        //   (3) Degree of Co-localization (PCC)
        //
        //Formula: Frequency Number * CoOccurrence+b1 * Degree of Co-Localization+b2
        //   b1 and b2 are baseline constants.  Here they are set to
        //   b1 = 0.5 and b2 = 0.2


        CoVals = cell(RowCounter,2);

        //titles
        CoVals(:,1) = CoLocalizations(:,1);
        CoVals{1,2} = 'CoVal';
        CoVals{1,3} = 'Frequency';
        CoVals{1,4} = 'CoDep}ence';
        CoVals{1,5} = 'CoLocalization';

        //CoVals
        for i = 2:RowCounter

            CoVals{i,3} = CoLocalizations{i,2}; //Frequency
            CoVals{i,4} = CoLocalizations{i,5}; //CoDep}ence
            CoVals{i,5} = CoLocalizations{i,8}; //CoLocalizations

           //             Frequency            * CoDep}ence
           CoVals{i,2} = (CoVals{i,3} * (CoVals{i,4} + 0.5) ...
        ... //   * CoLocalizations
               *(CoVals{i,5} + 0.2)^2)- 0.2;
        }

        //sort
        Header = CoVals(1,:);
        Body = sortrows(CoVals(2:RowCounter,:),2);

        CoVals = vertcat(Header,flipud(Body));

        }

        */
        return new CoRes();
    }

}
