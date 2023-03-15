    @staticmethod
    def statistics():
        dataset = []
        datasets = []
        subsets = []
        for item in os.listdir('./Peak Images/'):
            subsets = []
            if os.path.isdir('./Peak Images/' + item):
                for items in os.listdir('./Peak Images/' + item):
                    if items.endswith('.csv'):
                        with open('./Peak Images/'+ item + '/' + items, 'r') as file:
                            reader = csv.reader(file)
                            for row in reader:
                                if row[0] != "Mass/Charge Ratio":
                                    subsets.append([float(row[0]), float(row[1])])
                            for row in subsets:
                                datasets.append([row[0], row[1], item.split("_")[len(item.split("_"))-1] + " " + item.split("_")[len(item.split("_"))-6]])
                                dataset.append([row[0], row[1], item.split("_")[len(item.split("_"))-1] + " " + item.split("_")[len(item.split("_"))-6]])
        anovas = []
        kruskals = []
        mannwhit = []
        peaks = {}
        test = pd.DataFrame(datasets)
        test.columns =['Masses', 'Intensities', 'Target']
        lisss = datasets
        targerrr = test.loc[:,['Target']].values
        targerrrtwo = list(np.unique(targerrr))
        listofstuff = ['Masses']
        for val in targerrrtwo:
            listofstuff.append(str(val))
        testflipped = pd.DataFrame(columns=listofstuff)
        for index, row in test.iterrows():
            for columnName, columnData in testflipped.iteritems():
                if row['Target'] == columnName:
                    mass = row["Masses"]
                    intensity = row["Intensities"]
                    colnum = listofstuff.index(columnName)
                    rowid = len(testflipped.index)
                    filler = [row["Masses"]]
                    for num in range(0, len(targerrrtwo)):
                        filler.append(0)
                    testflipped.loc[len(testflipped.index)] = filler
                    testflipped.iloc[rowid][columnName] = intensity
        for index, row in testflipped.iterrows():
            for indextwo, rowtwo in testflipped.iterrows():
                if (row["Masses"] < rowtwo["Masses"] + .5) and (row["Masses"] > rowtwo["Masses"] - .5):
                    for col in row.iteritems():
                        if col[0] != "Masses":
                            if col[1] == 0 and testflipped.iloc[indextwo][col[0]] != 0:
                                testflipped.iloc[index][col[0]] = testflipped.iloc[indextwo][col[0]]
                                testflipped.iloc[index]["Masses"] = testflipped.iloc[indextwo]["Masses"]
        testflipped = testflipped.drop_duplicates()
        for row in subsets:
            num = round(row[0])
            if str(num) not in peaks:
                peaks[str(num)] = []
                for rows in subsets:
                    if((rows[0] > row[0] - 1) and (rows[0] < row[0] + 1)):
                        peaks[str(num)].append(rows[1])
                if(len(peaks[str(num)]) > 1):
                    #anovas.append(f_oneway(peaks[str(num)]))
                    #kruskals.append(kruskal(peaks[str(num)]))
                    #mannwhit.append(mannwhitneyu(peaks[str(num)], method="auto"))
                    pass
        peaklist = list(peaks.keys())
        """
        anovadata = {}
        anovadata['Peaks'] = peaklist
        anovadata['Anova'] = anovas
        kruskalsdata = {}
        kruskalsdata['Peaks'] = peaklist
        kruskalsdata['Kruskal Wallis'] = kruskals
        mannwhitdata = {}
        mannwhitdata['Peaks'] = peaklist
        mannwhitdata['Mann Whitney'] = mannwhit
        anovasdf = pd.DataFrame.from_dict(anovadata)
        kruskalsdf = pd.DataFrame.from_dict(kruskalsdata)
        mannwhitdf = pd.DataFrame.from_dict(mannwhitdata)
        """
        pcas = PCA(n_components=len(targerrrtwo)+1)
        features = listofstuff
        totaldatass = testflipped
        x = totaldatass.loc[:, features].values
        targerrr = listofstuff
        totaldatass.to_csv('Final Data Before PCA.csv', sep = '\t')
        x = StandardScaler(with_mean=False).fit_transform(x)
        principalComponents = pcas.fit_transform(x)
        pcas.fit(x)
        principalDf = pd.DataFrame(data = principalComponents, columns = targerrr)
        fig = plt.figure(figsize = (8,8))
        figtwo = plt.figure(figsize = (8,8))
        figthree = plt.figure(figsize = (8,8))
        elbowplot = figthree.add_subplot(1,1,1)
        scoreplot = figtwo.add_subplot(1,1,1)
        scoreplot.set_xlabel('Principal Component 1', fontsize = 15)
        scoreplot.set_ylabel('Principal Component 2', fontsize = 15)
        scoreplot.set_title('2 Component PCA', fontsize = 20)
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('Principal Component 1', fontsize = 15)
        ax.set_ylabel('Principal Component 2', fontsize = 15)
        ax.set_title('2 Component PCA', fontsize = 20)
        dictstuff = {}
        listofpoints = []
        targetcounter = []
        for target in targerrr:
            listofpointers = []
            if target != "Masses":
                for index, row in principalDf.iterrows():
                    for col in row.iteritems():
                        if col[0] == target:
                            listofpointers.append([row["Masses"], col[1]])
                x = []
                y = []
                xval = 0
                yval = 0
                for row in listofpointers:
                    x.append(row[0])
                    y.append(row[1])
                for item in x:
                    xval += item
                for item in y:
                    yval += item
                xval = xval/len(x)
                yval = yval/len(y)
                listofpoints.append([xval, yval])
                targetcounter.append([[xval, yval], target])
        import matplotlib.colors as mcolors
        colorlist = list(mcolors.BASE_COLORS)
        colornum = 6
        colors = []
        for num in range(0, colornum):
            colors.append(colorlist[num])
        """for target, colorval in zip(targerrr, colors):
            x = []
            y = []
            if target != "Masses":
                for coord in dictstuff[target]:
                    x.append(coord[0])
                    y.append(coord[1])
                scoreplot.scatter(x,y,color = colorval)
        """
        kmeans = KMeans(n_clusters=2)
        listofpointss = pd.DataFrame(listofpoints)
        q = kmeans.fit_predict(listofpointss)
        #filter rows of original data
        filtered_label0 = listofpointss[q == 1]
        filtered_label1 = listofpointss[q == 0]
        filtered_label0.columns = ["x", "y"]
        filtered_label1.columns = ["x", "y"]
        targlist1 = []
        targlist2 = []
        for index,row in filtered_label0.iterrows():
            mass = row[0]
            intensity = row[1]
            for rowtwo in targetcounter:
                    if rowtwo[1] != "Masses":
                        if rowtwo[0][0] == mass and rowtwo[0][1] == intensity:
                            targlist1.append(rowtwo[1])
        for index,row in filtered_label1.iterrows():
            mass = row[0]
            intensity = row[1]
            for rowtwo in targetcounter:
                    if rowtwo[1] != "Masses":
                        if rowtwo[0][0] == mass and rowtwo[0][1] == intensity:
                            targlist2.append(rowtwo[1])
        filtered_label1["Target"] = targlist2
        filtered_label0["Target"] = targlist1
        #plotting the results
        plottt = plt.figure(figsize = (8,8))
        scorer = plottt.add_subplot(1,1,1)
        scorer.set_xlabel('Principal Component 1', fontsize = 15)
        scorer.set_ylabel('Principal Component 2', fontsize = 15)
        scorer.set_title('2 Component PCA', fontsize = 20)
        for index,row in filtered_label0.iterrows():
            scorer.scatter(row["x"], row["y"], color = "red")
            scorer.annotate(row["Target"], xy = (row["x"], row["y"]))
        for index,row in filtered_label1.iterrows():
            scorer.scatter(row["x"], row["y"], color = "blue")
            scorer.annotate(row["Target"], xy = (row["x"], row["y"]))
        plottt.savefig("test.png")
        """
        targets = list((totaldatass.loc[:,['Target']].values).tolist())
        targs = []
        for targe in targets:
            if targe[0] not in targs:
                targs.append(targe[0])
        colornum = len(targs)
        colors = []
        finalDf = pd.concat([principalDf, totaldatass[['Target']]], axis = 1)
        import matplotlib.colors as mcolors
        colorlist = list(mcolors.BASE_COLORS)
        points = []
        targsss = []
        finalDf = principalDf
        for targer in targerrr:
            x = []
            y = []
            targ = ''
            for i in range(len(finalDf)):
                if finalDf.loc[i, 'Target'] == targer:
                    x.append(finalDf.loc[i, 'principal component 1'])
                    y.append(finalDf.loc[i, 'principal component 2'])
                    targ = targer
            pc1 = 0
            pc2 = 0
            for num in x:
                pc1 += num
            for num in y:
                pc2 += num
            pc1 = pc1/len(x)
            pc2 = pc2/len(y)
            points.append([pc1, pc2])
            targsss.append(targ)
        inertias = []
        newtargs = list(np.unique(targerrr))
        print(newtargs)
        q = len(newtargs) + 1
        for i in range(1,q):
            kmeans = KMeans(n_clusters=i)
            kmeans.fit(points)
            inertias.append(kmeans.inertia_)
        elbowplot.plot(range(1,q), inertias, marker='o')
        elbowplot.set_title('Elbow method')
        elbowplot.set_xlabel('Number of clusters')
        elbowplot.set_ylabel('Inertia')
        figthree.savefig('Elbow Grouping.png')
        pointss = pd.DataFrame(points)
        print(pointss)
        pointss.columns =['Masses', 'Intensities']
        kmeans = KMeans(n_clusters=2)
        q = kmeans.fit_predict(pointss)
        #filter rows of original data
        filtered_label0 = pointss[q == 1]
        filtered_label1 = pointss[q == 0]
        #plotting the results
        plottt = plt.figure(figsize = (8,8))
        scorer = plottt.add_subplot(1,1,1)
        scorer.scatter(filtered_label0['Masses'], filtered_label0['Intensities'], color = "red")
        scorer.scatter(filtered_label1['Masses'], filtered_label1['Intensities'], color = "blue")
        plottt.savefig('test.png')
        x = []
        y = []
        kmeans.fit(points)
        for point in points:
            x.append(point[0])
            y.append(point[1])
        colorss = kmeans.labels_
        for point in points:
            if colorss[points.index(point)] == 0:
                scoreplot.scatter(point[0], point[1], c = 'b', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
            if colorss[points.index(point)] == 1:
                scoreplot.scatter(point[0], point[1], c = 'g', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
            if colorss[points.index(point)] == 2:
                scoreplot.scatter(point[0], point[1], c = 'r', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
        for target, color in zip(targs,colors):
            indicesToKeep = finalDf['Target'] == target
            ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                    , finalDf.loc[indicesToKeep, 'principal component 2']
                    , c = color
                    , s = 50)
        ax.legend(['Delta', 'Omicron'])
        ax.grid()
        """
        figtwo.savefig("PCA Score Plot.png")
        fig.savefig('PCA Loading Plot.png')