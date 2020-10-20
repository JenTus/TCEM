#include "allocator.h"
#include "anyoption.h"
#include "iheap.h"
#include <ctime>
#include <numeric>
#include <iomanip>

namespace _Cide{
    
    allocator::allocator(AnyOption* opt1) {
        
        opt = opt1;
        probGraphFile = opt->getValue("x");
        compareFile = opt->getValue("p"); // it can be empty

        //delim = " \t";
        delim = " ";

        readGraphNodes(); // assign n and m

        sfmt_init_gen_rand(&sfmtSeed , 95082);

        if (compareFile == NULL) {
            cout << "Red k_r and k_b" << endl;
            k_r = strToInt(opt->getValue("r"));
            k_b = strToInt(opt->getValue("b"));
            cout << "kr is " << k_r << " " << "kb is " << k_b << endl;
        } else{
            std::string cF(compareFile);
            int beginIdx = cF.rfind('/');
            compareFileShort = cF.substr(beginIdx + 1);
            for (int j = 0; j < 4; ++j) {
                compareFileShort.pop_back();
            }
            readCompareNodes();
        }

        tao = ceil((double)k_b/(double)k_r);

//        P = Permutation(n, k_r*(tao + 1)) / (double) (pow(fact(tao), k_r) * fact(k_r));

        nrItems = 2;
        epsilon = strToDouble(opt->getValue("epsilon"));
        ell = strToDouble(opt->getValue("ell"));
        outFolderName = opt->getValue("outputFolder");
    
        nrPairs = n * nrItems; // multiple items

        nodeDegree.resize(n);
        
        for(int i = 0; i < n; i++)
            graphT.emplace_back(vector<int>());

        prevSize = 0; // when sample the first node
        //for(int i = 0; i < nrPairs; i++)
        //    hyperG.emplace_back(std::vector<int>());

        // create new pairs-related functions
        hyper_degreepairs = std::vector<double>(n*n, 0);
        for (int i = 0; i < n*n; ++i) {
           hyperGpairs.emplace_back(std::vector<int>());
        }

        rcList = new itemGraphList();
        //item-specific probT is aligned with graphT for each item
        for(int i = 0; i < nrItems; i++) {
            //_Cide::itemGraph *ig = new _Cide::itemGraph(n);
            auto *ig = new _Cide::itemGraph(n);
            for (int j = 0; j < n; j++)
                ig->probT.emplace_back(std::vector< double>());
            rcList->push_back(ig);
        }

        //nodeDegree = std::vector<double>(n,0.0);
        
        cout << "nr items " << nrItems << endl;
        //cout << "assignment size " << k << endl;

        readTICGraph();

        cout << "successfully read everything " << endl;

        arrangeOutputFiles();
        cout << "starting assignment for tcem" << endl;
        tdem();
        
// comment out below to run the baselines
        degreeVersionOne(theta);
        degreeVersionTwo(theta);
        mni(theta);
        if(compareFile!=NULL){
            compareGivenNodes(theta);
        }
    }


    double allocator::tdem() {
        
        cout << "computation for lower bounding OPT started " << endl;

        clock_t common_begin = clock();
        double lb = lowerBoundOPT();
        cout << "lower bound identified: " << lb << endl;
        
        theta = (2 + 2/3.0 * epsilon) * n * (ell * log(n) + log(2) +  logP(n,k_r,tao)) / (epsilon * epsilon * lb);
        //theta = (2 + 2/3.0 * epsilon) * n * (ell * log(n) + log(2) +  logcnk(n*(n-1), k_r*k_b)) / (epsilon * epsilon * lb);
        
        cout << "final sample size " << theta << endl;
        generateRCSets(theta);
        
        clock_t common_end = clock();
        duration_common = double(common_end - common_begin) / CLOCKS_PER_SEC;

        double greedySolution = n * rcGreedy(theta, true);

        return greedySolution; 
    }

    double allocator::lowerBoundOPT() {

        double epsilon_1 = epsilon * 2.0; //

        for (int x = 1; x < log2(n); x++) {

            int64 theta_x = (2+2.0/3.0 * epsilon_1)* (ell * log(n) + logP(n,k_r,tao) + log(log2(n))) * pow(2.0,x) / (epsilon_1 * epsilon_1);
            //int64 theta_x =  (2+2.0/3.0 * epsilon_1)* (ell * log(n) + logcnk(n*(n-1), k_r*k_b) + log(log2(n))) * pow(2.0,x) / (epsilon_1 * epsilon_1);
            
//            std::cout << "theta_x " << theta_x << std::endl;
            generateRCSets(theta_x);
            
            double ept = rcGreedy(theta_x, false);

//            cout << "ept is: " << ept << endl;

            if (ept > ((1+epsilon_1) / pow(2.0, x))) {
                double lowerBound = ept * n / (1 + epsilon_1);
                // kontrol:
//                std:cout << "x and theta_x " << x << " " << theta_x << std::endl;
                return lowerBound;
            }
        }
        cout << "returning naive lower bound  " << endl;
        double naive = 1.0;
        return naive;
    }

    void allocator::generateRCSets(int64 newSize) {
        
        if(newSize < prevSize){
            theta = prevSize; //
            return;
        }

        // sample target nodes
        for (int i = prevSize; i < newSize; i++) {
            int randTarget = sfmt_genrand_uint32(&sfmtSeed) % n;
            targetNodes.push_back(randTarget);
        }
        
        // expand coordinated RR sets samples of items
        for(int itemID = 0; itemID < nrItems; itemID++) {
            rcList->at(itemID)->generateRRSample(targetNodes, prevSize, newSize);
        }

        for (int rcID = prevSize; rcID < newSize; ++rcID) {
            int tprsize = rcList->at(0)->hyperGT[rcID].size();
            for (int j = 0; j < tprsize; ++j) {
                int y = rcList->at(0)->hyperGT[rcID][j];
                //cout << y << " ";
            }
            int tpbsize = rcList->at(1)->hyperGT[rcID].size();
            for (int j = 0; j < tpbsize; ++j) {
                int y = rcList->at(1)->hyperGT[rcID][j];
                //cout << y << " ";
            }
            //cout << endl;
        }

        // rcList -> itemID -> hyperGT[rcID]: store a list of the nodes.
        double tmpsize = createhyperGTPairs(newSize);
//        cout << "current size is " << newSize <<  " size of hyPerGTPairs is " << tmpsize << endl;

        for (int rcID = prevSize; rcID < newSize; ++rcID) {
            int tpsize = hyperGTpairs[rcID].size();
            for (int i = 0; i < tpsize; ++i) {
                int x = hyperGTpairs[rcID][i]; // because x is too large.
                hyperGpairs[x].push_back(rcID);
                hyper_degreepairs[x] += 1;
            }
            //int tprsize = rcList->at(0)->hyperGT[rcID].size();
            //for (int j = 0; j < tprsize; ++j) {
            //    int y = rcList->at(0)->hyperGT[rcID][j];
            //    //cout << y << " ";
            //    hyperG[y].push_back(rcID); // careful about the repetition
            //}
            //int tpbsize = rcList->at(1)->hyperGT[rcID].size();
            //for (int j = 0; j < tpbsize; ++j) {
            //    int y = rcList->at(1)->hyperGT[rcID][j];
            //    //cout << y << " ";
            //    hyperG[y+n].push_back(rcID); // careful about the repeatation
            //}
            //cout << endl;
        }

        prevSize = newSize;
        
    }

    double allocator::createhyperGTPairs(int64 newSize) {
        // use to create reverse rr pairs
        //update hyperG and hyper_degree
        //cout << "newSize is " << newSize << endl;
        for (int z = prevSize; z < newSize; ++z) {
            hyperGTpairs.emplace_back(std::vector<int>()); // init
        }

        for (int rcID = prevSize; rcID < newSize; ++rcID) {
            for (int i = 0; i < rcList->at(0)->hyperGT[rcID].size(); ++i) {
                for (int j = 0; j < rcList->at(1)->hyperGT[rcID].size(); ++j) {
                    //cout << "r: " << rcList->at(0)->hyperGT[rcID].size() << " b: " << rcList->at(1)->hyperGT[rcID].size()<< endl;
                   int t = rcList->at(0)->hyperGT[rcID][i]*n + rcList->at(1)->hyperGT[rcID][j];
                   hyperGTpairs[rcID].push_back(t);
                }
            }
        }

        return hyperGTpairs.size();
    }


    double allocator::rcGreedy(int64 rcSampleSize, bool extraResults) {
        // For initialization the k_r and k_b, to delete later.

        double totalExpScore = 0;
        std::vector<double> ExpScore(n*n, 0); // use to store hyper_degreepairs

        int bestPairID = -1, bestRedNodeID = -1, bestBlueNodeID = -1;
        std::vector<int> isValidRed = std::vector<int>(n,tao);
        std::vector<int> isValidBlue = std::vector<int>(n,1);
        std::vector<bool> isCovered = std::vector<bool>(rcSampleSize, false);
        std::set<int> rednodes;
        std::set<int> bluenodes;

        int flag_red = 0;

        clock_t begin = clock();

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (i != j){
                    ExpScore[i*n + j] = hyper_degreepairs[i*n + j];
                }
            }
        }

        seedSet.clear();
        seedScores.clear(); // seedScores should be cleared as well.

        while ((int)seedSet.size() < k_b) {

            double tempsum = accumulate(ExpScore.begin(), ExpScore.end(), 0.01);
            if(tempsum < 1){
                break;
            }
            // here I use the vector to store the results, but can use the priority queue.
            //bestPairID = std::max_element(ExpScore.begin(), ExpScore.end()) - ExpScore.begin();
            bestPairID = std::distance(ExpScore.begin(), std::max_element(ExpScore.begin(), ExpScore.end()));
            //cout << "bestPairID is " << bestPairID   <<endl;
            bestRedNodeID = (int) bestPairID / n;
            rednodes.insert(bestRedNodeID);
            bestBlueNodeID = (int) bestPairID % n;
            bluenodes.insert(bestBlueNodeID);

            // if can not select, then just pop out
            if (isValidRed[bestRedNodeID] == 0 || isValidBlue[bestBlueNodeID] == 0) {
                ExpScore[bestPairID] = 0;
                continue;
            }

            // push
            seedSet.push_back(bestPairID); // store
            //cout << "ExpScore is " << ExpScore[bestPairID] << "rc sample size is " <<rcSampleSize << "n is " << n << endl;
            seedScores.push_back(ExpScore[bestPairID]); // store score, dividing later

            // make update
            isValidRed[bestRedNodeID] -= 1;
            isValidRed[bestBlueNodeID] = 0;
            isValidBlue[bestBlueNodeID] = 0;
            isValidBlue[bestRedNodeID] = 0;

            if(isValidRed[bestRedNodeID] <= 0){
                for (int i = 0; i < n; ++i) {
                   ExpScore[bestRedNodeID*n + i] = 0;
                }
            }

            for (int i = 0; i < n; ++i) {
               ExpScore[i*n + bestBlueNodeID]  = 0;
               ExpScore[i*n + bestRedNodeID]  = 0;
               ExpScore[bestBlueNodeID*n + i] = 0;
            }

            // if red nodes reaches k_r
            if(rednodes.size() == k_r && flag_red == 0){
                for (int i = 0; i < n; ++i) {
                    if(rednodes.count(i) == 0){
                        isValidRed[i] = 0;
                        for (int j = 0; j < n; ++j) {
                           ExpScore[i*n + j] = 0;
                        }
                    }
                }
                flag_red = 1;
            }

//            cout << "sum of Expscore is " << accumulate(ExpScore.begin(), ExpScore.end(), 0) << endl;

            // update the ExpScore of the pairs whose marginal gain change due to the best pair selection
            int tmpsize = hyperGpairs[bestPairID].size();
            for (int j = 0; j < tmpsize; j++) {
                int rcID = hyperGpairs[bestPairID][j]; // get the rcID
                if(isCovered[rcID] == false){ // ensure not to repeatedly --
                    int ttmpsize = hyperGTpairs[rcID].size();
                    for (int z = 0; z < ttmpsize; z++) {
                        int pairID = hyperGTpairs[rcID][z];
                        if (ExpScore[pairID] != 0) {
                            ExpScore[pairID]--;
                        }
                    }
                    isCovered[rcID] = true;
                }
           }
        } // end of while k loop

        for (int i1 = 0; i1 < seedScores.size(); ++i1) {
           totalExpScore += seedScores[i1];
        }

        int covercount = 0;
        for (int l = 0; l < isCovered.size(); ++l) {
            if(isCovered[l])
                covercount++;
        }


        ///////////// produce the additional results here
        ////////////  Remember should write the f function results not the g function result.
        if(extraResults) {

        cout << "totalExpScore is: " << totalExpScore << " cover is " << covercount << " sample size is " << rcSampleSize << endl;

            clock_t end = clock();

            //cout << "total tdem time taken (in seconds) " << totalDuration << " total memory (in mb) " << totalMemory << endl;
            int redbluecount = 0;

            vector<int> s1(rednodes.begin(), rednodes.end());
            vector<int> s2(bluenodes.begin(), bluenodes.end());
            sort(s1.begin(), s1.end());
            sort(s2.begin(), s2.end());


            for (int i = 0; i < rcSampleSize; ++i) {
                vector<int> RedList = rcList->at(0)->hyperGT[i];
                vector<int> BlueList = rcList->at(1)->hyperGT[i];

                sort(RedList.begin(), RedList.end());
                sort(BlueList.begin(), BlueList.end());

                if(intersection_size(s1, RedList) >= 1){
                    if(intersection_size(s2, BlueList) >= 1){
                        redbluecount +=1;
                    }
                }
            }

            totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
            totalMemory = getCurrentMemoryUsage(); // in MB

            float totalscore = (float) redbluecount / (float) rcSampleSize;

            cout << "Total score of Greedy is " << totalscore*n << "  g result is " << totalExpScore / (float) rcSampleSize *n  << endl;

            std::ostringstream streamObj1;
            streamObj1 << std::fixed << std::setprecision(2)<<totalscore*n << "(" << totalExpScore/(float) rcSampleSize*n << ")";

            writeInMasterOutputFile("Greedy", totalDuration,  totalMemory, streamObj1.str(), rednodes, bluenodes);

        } // end of extra results
        
        return totalExpScore / (double)rcSampleSize;
    }


    double allocator::degreeVersionOne(int64 rcSampleSize) {

        clock_t begin = clock();

        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond> heap;

        float totalscore = 0;
        int bestNodeID;
//
//
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double> pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }

        seedSet.clear();
        seedScores.clear();
        set<int> rednodes;
        set<int> bluenodes;
        int testt;


        while ((int) seedSet.size() < k_r) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            seedSet.push_back(bestNodeID);
            rednodes.insert(bestNodeID);
        }


        while ((int) seedSet.size() < k_r + k_b) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
            seedSet.push_back(bestNodeID);
            bluenodes.insert(bestNodeID);
        }

        int redbluecount = 0;

        vector<int> s1(rednodes.begin(), rednodes.end());
        vector<int> s2(bluenodes.begin(), bluenodes.end());
        sort(s1.begin(), s1.end());
        sort(s2.begin(), s2.end());


        for (int i = 0; i < rcSampleSize; ++i) {
            vector<int> RedList = rcList->at(0)->hyperGT[i];
            vector<int> BlueList = rcList->at(1)->hyperGT[i];

            sort(RedList.begin(), RedList.end());
            sort(BlueList.begin(), BlueList.end());

            if(intersection_size(s1, RedList) >= 1){
                if(intersection_size(s2, BlueList) >= 1){
                    redbluecount +=1;
                }
            }
        }

        //cout << " redBluecount is " << redbluecount << endl;

        totalscore = (float) redbluecount/ (float) rcSampleSize;

        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB

        cout << "Total score of degree one is " << totalscore*n << endl;

        std::ostringstream streamObj1;
        streamObj1 << std::fixed << std::setprecision(2)<<totalscore*n;

        writeInMasterOutputFile("DegreeVersionOne", totalDuration,  totalMemory, streamObj1.str(), rednodes, bluenodes);

        return totalscore;
    }


    double allocator::degreeVersionTwo(int64 rcSampleSize) {


        clock_t begin = clock();
        priority_queue<pair<int, double>, vector<pair<int, double>>, CompareBySecond> heap;

        float totalscore = 0;
        int bestNodeID;
//
        for (int i = 0; i < n; i++) {
            double val = (double) nodeDegree[i];
            std::pair<int, double> pairVal(std::make_pair(i, val));
            heap.push(pairVal);
        }

        std::set<int> rednodes, bluenodes;
        seedSet.clear();
        seedScores.clear();
        //std::set<int> rrRed;
        //std::set<int> rrBlue;
        //int testt;

        //rrRed.reserve(n*k_r);
        //rrBlue.reserve(n*k_b);

        while ((int) seedSet.size() < k_r + k_r) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
           // for (int j = 0; j < hyperG[bestNodeID].size(); ++j) {
           //     int id = hyperG[bestNodeID][j];
           //     rrRed.insert(id);
           // }

            seedSet.push_back(bestNodeID);
            rednodes.insert(bestNodeID);

            pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
           // for (int j = 0; j < hyperG[bestNodeID + n].size(); ++j) {
           //     int id = hyperG[bestNodeID + n][j];
           //     rrBlue.insert(id);
           // }

            seedSet.push_back(bestNodeID);
            bluenodes.insert(bestNodeID);
        }

        while ((int) seedSet.size() < k_r + k_b) {
            pair<int, double> pairVal = heap.top();
            heap.pop();
            bestNodeID = pairVal.first;
           // for (int l = 0; l < hyperG[bestNodeID + n].size(); ++l) {
           //     int id = hyperG[bestNodeID + n][l];
           //     rrBlue.insert(id);
           // }
            seedSet.push_back(bestNodeID);
            bluenodes.insert(bestNodeID);
        }

        /*

        vector<int> rrRedV(rrRed.begin(), rrRed.end());
        vector<int> rrBlueV(rrBlue.begin(), rrBlue.end());

        std::sort(rrRedV.begin(), rrRedV.end());
        std::sort(rrBlueV.begin(), rrBlueV.end());

        testt = intersection_size(rrRedV, rrBlueV);
        cout << "intersection is " << testt;
        */

        int redbluecount = 0;

        vector<int> s1(rednodes.begin(), rednodes.end());
        vector<int> s2(bluenodes.begin(), bluenodes.end());
        sort(s1.begin(), s1.end());
        sort(s2.begin(), s2.end());


        for (int i = 0; i < rcSampleSize; ++i) {
            vector<int> RedList = rcList->at(0)->hyperGT[i];
            vector<int> BlueList = rcList->at(1)->hyperGT[i];

            sort(RedList.begin(), RedList.end());
            sort(BlueList.begin(), BlueList.end());

            if(intersection_size(s1, RedList) >= 1){
                if(intersection_size(s2, BlueList) >= 1){
                    redbluecount +=1;
                }
            }
        }


        //cout << " redBluecount is " << redbluecount << endl;

        totalscore = (float) redbluecount / (float) rcSampleSize;

        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB


        cout << "Total score of degree two is " << totalscore*n << endl;
        std::ostringstream streamObj1;
        streamObj1 << std::fixed << std::setprecision(2)<<totalscore*n;

        writeInMasterOutputFile("DegreeVersionTwo", totalDuration,  totalMemory, streamObj1.str(), rednodes, bluenodes);
        return totalscore;
    }



    double allocator::compareGivenNodes(int64 rcSampleSize) {


        clock_t begin = clock();

        float totalscore = 0;
        int bestNodeID;
        seedSet.clear();
        seedScores.clear();


        int redbluecount = 0;

        sort(comparered.begin(), comparered.end());
        sort(compareblue.begin(), compareblue.end());


        for (int i = 0; i < rcSampleSize; ++i) {
            vector<int> RedList = rcList->at(0)->hyperGT[i];
            vector<int> BlueList = rcList->at(1)->hyperGT[i];

            sort(RedList.begin(), RedList.end());
            sort(BlueList.begin(), BlueList.end());

            if(intersection_size(comparered, RedList) >= 1){
                if(intersection_size(compareblue, BlueList) >= 1){
                    redbluecount +=1;
                }
            }
        }

        //cout << " redBluecount is " << redbluecount << endl;
        totalscore = (float) redbluecount / (float) rcSampleSize;
        // end of extra results
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB


        cout << "Total score of BalanceExposure is " << totalscore*n << endl;
        std::ostringstream streamObj1;
        streamObj1 << std::fixed << std::setprecision(2)<<totalscore*n;

        set<int> setcomparered(comparered.begin(), comparered.end());
        set<int> setcompareblue(compareblue.begin(), compareblue.end());

        writeInMasterOutputFile("BalanceExposure", totalDuration,  totalMemory, streamObj1.str(), setcomparered, setcompareblue);
        return totalscore;
    }

    double allocator::mni(int64 rcSampleSize){
        clock_t begin = clock();

        std::vector<std::vector<int>> graph; // the graph is created from graphT
        std::vector<std::vector<int>> graphPairs; // for each pair, store common nodes
        std::vector<int> graphPairsDegree(n*n, 0); //
        std::vector<std::vector<int>> graphTPairs; // transpose of the graph, for node i, store the pairs that cover i
        graphPairs.reserve(n*n);
        graphTPairs.reserve(n);
        graph.reserve(n);
        for (int i = 0; i < n; ++i) {
           graph.emplace_back(vector<int>());
           graphTPairs.emplace_back(vector<int>());
        }

        for (int k1 = 0; k1 < n*n; ++k1) {
            graphPairs.emplace_back(vector<int>());
        }

        //create graph
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < graphT[j].size(); ++i) {
               graph[graphT[j][i]].push_back(j);
            }
        }
        //for test

        for (int i1 = 0; i1 < n; ++i1) {
           graph[i1].push_back(i1);
        }

        //create graphPairs
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                if(i == k)
                    continue;
                std::vector<int> tmpv;
                tmpv = intersection(graph[k], graph[i]);
                graphPairs[k*n + i].insert(graphPairs[k*n +i].begin(), tmpv.begin(), tmpv.end());
                graphPairsDegree[k*n + i] = graphPairs[k*n + i].size(); //degree

                for (int j = 0; j < graphPairs[k*n + i].size(); ++j) {
                    graphTPairs[graphPairs[k*n + i][j]].push_back(k*n + i);
                }
            }
        }

        int bestPairID, bestRedNodeID, bestBlueNodeID;
        std::vector<int> isValidRed = std::vector<int>(n,tao);
        std::vector<int> isValidBlue = std::vector<int>(n,1);
        std::vector<bool> isCovered = std::vector<bool>(n, false);
        seedSet.clear();
        seedScores.clear();
        double totalExpScore = 0;
        std::set<int> rednodes;
        std::set<int> bluenodes;
        int flag_red = 0;

        while ((int)seedSet.size() < k_b) {
            double tempsum = accumulate(graphPairsDegree.begin(), graphPairsDegree.end(), 0.01);
            if(tempsum < 1){
                break;
            }

            bestPairID = std::distance(graphPairsDegree.begin(), std::max_element(graphPairsDegree.begin(), graphPairsDegree.end()));
            //cout << "bestPairID is " << bestPairID   <<endl;
            bestRedNodeID = (int) bestPairID / n;
            rednodes.insert(bestRedNodeID);
            bestBlueNodeID = (int) bestPairID % n;
            bluenodes.insert(bestBlueNodeID);

            //cout << "sum of pairwise common neighbors is " << accumulate(graphPairsDegree.begin(), graphPairsDegree.end(), 0) << endl;

            if (isValidRed[bestRedNodeID] == 0 || isValidBlue[bestBlueNodeID] == 0) {
                graphPairsDegree[bestPairID] = 0;
                continue;
            }

            seedSet.push_back(bestPairID); // store
            seedScores.push_back(graphPairsDegree[bestPairID]); // store score, dividing later

            // make update
            isValidRed[bestRedNodeID] -= 1;
            isValidRed[bestBlueNodeID] = 0;
            isValidBlue[bestBlueNodeID] = 0;
            isValidBlue[bestRedNodeID] = 0;

            if(isValidRed[bestRedNodeID] <= 0){
                for (int i = 0; i < n; ++i) {
                    graphPairsDegree[bestRedNodeID*n + i] = 0;
                }
            }

            for (int i = 0; i < n; ++i) {
                graphPairsDegree[i*n + bestBlueNodeID]  = 0;
                graphPairsDegree[i*n + bestRedNodeID]  = 0;
                graphPairsDegree[bestBlueNodeID*n + i] = 0;
            }

            // if red nodes reaches k_r
            if(rednodes.size() == k_r && flag_red == 0){
                for (int i = 0; i < n; ++i) {
                    if(rednodes.count(i) == 0){
                        isValidRed[i] = 0;
                        for (int j = 0; j < n; ++j) {
                            graphPairsDegree[i*n + j] = 0;
                        }
                    }
                }
                flag_red = 1;
            }

            // update the ExpScore of the pairs whose marginal gain change due to the best pair selection
            int tmpsize = graphPairs[bestPairID].size();
            for (int j = 0; j < tmpsize; j++) {
                int nodeid = graphPairs[bestPairID][j];
                if(isCovered[nodeid] == false){ // ensure not to repeatedly --
                    int ttmpsize = graphTPairs[nodeid].size();
                    for (int z = 0; z < ttmpsize; z++) {
                        int pairID = graphTPairs[nodeid][z];
                        if (graphPairsDegree[pairID] != 0) {
                            graphPairsDegree[pairID]--;
                        }
                    }
                    isCovered[nodeid] = true;
                }
            }
        } // end of while k loop

        for (int i1 = 0; i1 < seedScores.size(); ++i1) {
            totalExpScore += seedScores[i1];
        }

        int covercount = 0;
        for (int l = 0; l < isCovered.size(); ++l) {
            if(isCovered[l])
                covercount ++;
        }

        cout << "number of common covered neighbor is: " << totalExpScore << " cover is " << covercount << endl;

        // get the score on rr sets
        clock_t end = clock();
        totalDuration = double(end - begin) / CLOCKS_PER_SEC + duration_common; // in seconds
        totalMemory = getCurrentMemoryUsage(); // in MB


        int redbluecount = 0;

        vector<int> s1(rednodes.begin(), rednodes.end());
        vector<int> s2(bluenodes.begin(), bluenodes.end());
        sort(s1.begin(), s1.end());
        sort(s2.begin(), s2.end());


        for (int i = 0; i < rcSampleSize; ++i) {
            vector<int> RedList = rcList->at(0)->hyperGT[i];
            vector<int> BlueList = rcList->at(1)->hyperGT[i];

            sort(RedList.begin(), RedList.end());
            sort(BlueList.begin(), BlueList.end());

            if(intersection_size(s1, RedList) >= 1){
                if(intersection_size(s2, BlueList) >= 1){
                    redbluecount +=1;
                }
            }
        }

        float totalscore = (float) redbluecount / (float) rcSampleSize;

//
        cout << "Total score of mni is " << totalscore*n << endl;

        std::ostringstream streamObj1;
        streamObj1 << std::fixed << std::setprecision(2)<<totalscore*n << "(" << totalExpScore << ")";

        writeInMasterOutputFile("MNI", totalDuration,  totalMemory, streamObj1.str(), rednodes, bluenodes);

        return totalscore*n;

    }


    void allocator::readGraphNodes() {

        //string probGraphFile = opt->getValue("probGraphFile");
        ifstream myfile (probGraphFile.c_str(), ios::in);

        //double *itemProbs; //item-specific influence probabilities

        int nrEdges = 0;
        int nrNodes = 0;

        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;

                std::vector<std::string> nodeids;
                split1(line, nodeids);
                 int u1 = strToInt(nodeids[0]);
                 int u2 = strToInt(nodeids[1]);
            
                if (u1 == u2)
                    continue;

                nrEdges++;
                // for control

                if(u1>nrNodes){
                    nrNodes = u1;
                }
                if(u2>nrNodes){
                    nrNodes = u2;
                }

            }

            myfile.close();
        }

        n = nrNodes +1;
        m = nrEdges;
        cout << "total number of nodes " << n << endl;
        cout << "total number of edges " << nrEdges << endl;
    }

    void allocator::readCompareNodes() {

        //string probGraphFile = opt->getValue("probGraphFile");
        ifstream myfile (compareFile, ios::in);

        //double *itemProbs; //item-specific influence probabilities

        if (myfile.is_open()) {
            string line1;
            string line2;
            getline(myfile, line1);
            getline(myfile, line2);

            vector<string> nodeids1;
            vector<string> nodeids2;

            split1(line1, nodeids1);
            split1(line2, nodeids2);

            if(nodeids1.size() > nodeids2.size()){
                for (int i = 2; i < nodeids1.size(); ++i) {
                    unsigned int tt = strToInt(nodeids1[i]);
                    compareblue.push_back(tt);
                }
                for (int i = 2; i < nodeids2.size(); ++i) {
                    unsigned int tt = strToInt(nodeids2[i]);
                    comparered.push_back(tt);
                }
            }
            else{
                for (int i = 2; i < nodeids1.size(); ++i) {
                    unsigned int tt = strToInt(nodeids1[i]);
                    comparered.push_back(tt);
                }
                for (int i = 2; i < nodeids2.size(); ++i) {
                    unsigned int tt = strToInt(nodeids2[i]);
                    compareblue.push_back(tt);
                }
            }

            myfile.close();
        }

        k_r = comparered.size();
        k_b = compareblue.size();
        cout << "total number of red nodes " << k_r << endl;
        cout << "total number of blue nodes " << k_b << endl;
    }


    void allocator::readTICGraph() {
        
        //string probGraphFile = opt->getValue("probGraphFile");
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        //double *itemProbs; //item-specific influence probabilities

        int nrEdges = 0;
        set<int> nodes; // for control
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;

                std::vector<std::string> nodeids;
                split1(line, nodeids);
                 int u1 = strToInt(nodeids[0]);
                 int u2 = strToInt(nodeids[1]);
                double u3 = strToDouble(nodeids[2]); // the probability
                double u4 = strToDouble(nodeids[3]);
                
                if (u1 == u2)
                    continue;

                nodeDegree[u1] = nodeDegree[u1] + 1.0;

                nrEdges++;
                
                graphT[u2].push_back(u1); //insert to the transposed graph
                
                // for control
                nodes.insert(u1);
                nodes.insert(u2);

                rcList->at(0)->probT[u2].push_back(u3); // insert the probabilty
                rcList->at(1)->probT[u2].push_back(u4);

            }


            myfile.close();
        }
        
        else
            cout << "Can't open input graph file " << probGraphFile << endl;
        
        cout << "graph import complete " << endl;
        
    }

    
    allocator::~allocator() {
        cout << "assignments complete! " << endl;
    }
    
    void allocator::arrangeOutputFiles() {
        
       string command = string("mkdir -p ") + outFolderName ;
        
       system(command.c_str());

        std::size_t pos = probGraphFile.find("/");
        std::size_t end = probGraphFile.find(".txt");

        string masterFileName = probGraphFile.substr(pos+1, end-pos-1) + "+" + to_string(k_r) + "+" + to_string(k_b) + ".txt";

        outMasterName = outFolderName + OS_SEP + masterFileName;
//        outMasterName = masterFileName;
        
        if(outMasterStream.is_open())
            outMasterStream.close();
        
        outMasterStream.open(outMasterName.c_str());
        
        outMasterStream << "n" << " " << "m" << " " << "redNodes" << " " << "blueNodes" << " "<< endl;
        outMasterStream << n << " " << m << " " << k_r << " " << k_b << " " << endl;
        outMasterStream << "Algorithm" << " " << "Time" << " " << "Memory(MB)" << " " << "score" << endl;
        
        if (!outMasterStream.is_open()) {
            cout << "Can't open file " << outMasterName  << " for writing" << endl;
            exit(1);
        }
        
        // memory icin
        command = string("mkdir -p temp") ;
        system(command.c_str());
    }
    
    //    void allocator::writeInMasterOutputFile(int nodeID, int itemID) {
    //        // seed-node item mgScore totScore runTime(sec) memory(mb)
    //        outMasterStream << nodeID << " " << itemID << endl;
    //    }
    
    
    void allocator::writeInMasterOutputFile(const string& algorithm, float duration, float memory, string score, set<int> reds, set<int> blues) {
        // seed-node item mgScore totScore runTime(sec) memory(mb)
        outMasterStream << std::fixed<<std::setprecision(2) << algorithm << " & " << duration << " & " <<  memory << " & " << score <<  endl;
        set<int>::iterator it;

        outMasterStream << "red: ";
        for (it=reds.begin(); it!=reds.end(); ++it)
            outMasterStream << *it << " ";
        outMasterStream << endl;

        outMasterStream << "blue: ";
        for (it=blues.begin(); it!= blues.end(); ++it) {
            outMasterStream << *it << " ";
        }
        outMasterStream << endl;
    }


}
