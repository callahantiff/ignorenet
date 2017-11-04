################################################################################################
# AnnotationProcessing.py
# Purpose: script processes the annotation results from PubAnnotation
# version 1.0.0
# date: 10.24.2017
################################################################################################


## import module/script dependencies
import QueryEndpoint
import json



def main():

    # read in query
    file = 'Queries/preeclampsia_10.25.17.txt'
    Query = open(file).read()

    # run query against endpoint
    url = 'http://pubannotation.org/projects/Preeclampsia/search'
    output = 'Queries/Results/preeclampsia_10.25.17'
    annotations = QueryEndpoint.RunQuery(Query, output, url)
    len(annotations['results']['bindings'])

    # process results
    result_dict = {}
    for result in annotations['results']['bindings']:
        print('\n')
        for var in annotations['head']['vars']:
            print(var, result[var]['value'].split('/')[-1])

            # label keys/values
            key = var
            value = result[var]['value'].split('/')[-1]

            if var in result_dict.keys():
                result_dict[key].append(value)

            else:
                result_dict[key] = [value]











if __name__ == "__main__":
    main()