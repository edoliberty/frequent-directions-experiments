import sys, json, numpy, traceback
from time import time as timer
from dataHandler import DataHandler
from matrixSketcher import MatrixSketcher
from matrixSketcherSVD import MatrixSketcherSVD
from matrixSampler import MatrixSampler
from matrixMultiplier import MatrixMultiplier
from matrixNothing import MatrixNothing
from matrixRandomSums import  MatrixRandomSums
from matrixRandomProjections import MatrixRandomProjections
#import xlwt
import csv

class DotPrinter():
    def __init__(self,n1=1,n2=100):
        self.i = 0
        self.n1 = n1
        self.n2 = n2
        self.bsn = False
    def last(self):
        self.i = 0
        return '' if self.bsn else '\n'
        
    def next(self):        
        self.i += 1
        s = ''
        if i%self.n1==0:
            s = '.'            
        if self.i >= self.n2:
            self.i=0
            s+='\n'
            self.bsn = True
        else:
            self.bsn = False    
        return s  
        
        
def listOuterProduct(single_item_lists):
    op_list = [[]]
    for single_item_list in single_item_lists:
        temp_op_list = []
        for op_item in op_list:
            for single_item in single_item_list:
                new_item = [x for x in op_item] + [single_item]
                temp_op_list.append(new_item)
        op_list = temp_op_list 
    #print 'op_list = ',op_list
    return op_list
    
def getCommandLineArguments(cammand_line_arguments):
    params = {}
    params['N'] = [10000]#range(10000,110000,10000)#[10000]
    params['round'] = [1]
    params['dimension'] = [1000]
    params['ell'] = range(10,220,10)
    params['signal_dimension'] = [50]
    params['signal_to_noise_ratio'] = [10]
    params['signal_singular_value_decay_factor'] = [1]
    params['signal_singular_value_decay_type'] = ['lin']#'exp','lin']
    
    for arg in cammand_line_arguments:
        arg_parts = arg.split('=')        
        if not len(arg_parts) == 2:
            continue
        (key,values_string) = arg_parts
        key = key.strip('-')
        values = [int(v) for v in values_string.strip('()[]{}').split(',')]
        params[key] = values
    
    params_items = sorted(params.items())
    params_names = [item[0] for item in params_items]
    params_lists = [item[1] for item in params_items]
    #print 'here'
    #print params_names
    #print params_lists
    return params_names, params_lists
    

                    
if __name__ == '__main__':    
    run_mode = not '-test' in sys.argv
    params_names, params_lists = getCommandLineArguments(sys.argv[1:])      
    
    scketcherClasses = [MatrixMultiplier, MatrixSketcher, MatrixSampler, MatrixRandomSums, MatrixRandomProjections, MatrixNothing]
    #scketcherClasses = [MatrixNothing,MatrixMultiplier]

    dp = DotPrinter(100,1000) 
    results_computed = 0
    
    outfile = open('projected_sketch_d_50.csv','w')
    outwriter = csv.writer(outfile)
    row_to_write = ["MatrixMultiplier", "MatrixSketcher", "MatrixSampler", "MatrixRandomSums", "MatrixRandomProjections", "MatrixNothing"]
    data_to_write = []

    for params_values in listOuterProduct(params_lists):
        data_to_write.append(row_to_write)
        try:
            # this is a hacky way to set all the locals to their values
            
            local_param_dict = dict(zip(params_names,params_values))
            locals().update(local_param_dict)
            #value = ('ell = ',ell,'signal_dimension = ',signal_dimension)
            #print 'val = ',str(value)
            #outfile.write(str(value))
            #outfile.write('\n')
            if ell > dimension/2:
                sys.stderr.write('INFO: Not computing parameters %s\n'%str(local_param_dict))
                continue
            else:
                sys.stderr.write('INFO: Computing parameters %s\n'%str(local_param_dict))
            
            # remove the brute force O(d^2) computation when the dimension is too large  
            dimension_too_large_for_brute_force = dimension > 3000
            if dimension_too_large_for_brute_force:
                scketcherClasses.remove(MatrixMultiplier)
                sys.stderr.write('INFO: Not including MatrixMultiplier since the dimension is %d\n'%dimension)
            
            # initializing sketches
            sketchers = [sketch_class(dimension,ell,ell/2) for sketch_class in scketcherClasses]
            #print 'after initializing sketchers'
            
            # initializing data handler
            sys.stderr.write('INFO: Initializing data handler\n')
            data_handler = DataHandler()
            data_handler.initBeforeMake(dimension,signal_dimension,signal_to_noise_ratio,\
                                        signal_singular_value_decay_factor,signal_singular_value_decay_type)
            #print 'after dataHandler'
            update_times = dict([(sketcher.class_name,0.0) for sketcher in sketchers])
            
            sys.stderr.write('INFO: Updating sketchers\n')
            last_id_dot = False
            #try:
            A = numpy.zeros((N,dimension))
            if run_mode:
                for i in xrange(N):
                    v = data_handler.makeRow()
                    A[i,:] = v
                    for (j,sketcher) in enumerate(sketchers):
                        t_start = timer()
                        sketcher.add(v)
                        t_end = timer()
                        update_times[sketcher.class_name] += t_end-t_start
                    
                    sys.stderr.write(dp.next())
                sys.stderr.write(dp.last())
            sys.stderr.write('INFO: Finished updating sketchers\n')
            
            
            # Getting the exact computation 
            if not dimension_too_large_for_brute_force:
                for sketcher in sketchers:
                    if sketcher.class_name == 'MatrixMultiplier':
                        ATA = sketcher.getATA()
            
            signal_row_space = data_handler.getSignalRowSpace()
            
            #print 'before sketchers for loop'
            row_to_write = []

            for sketcher in sketchers:
                total_res = 0
                for itr in range(5):
                    sys.stderr.write('INFO: Computing results for %s\n'%sketcher.class_name)
                    res = local_param_dict.copy()
                    res['sketcher_class_name'] = sketcher.class_name
                
                # getting the sketch + performance metrics
                    res['update_time'] = update_times[sketcher.class_name]
                    t_start = timer()
                    sketch = sketcher.get()
                    t_end = timer()
                    res['sketch_compute_time'] = t_end-t_start
                    res['total_time'] = res['sketch_compute_time'] + res['update_time']  
                
            
                # quality metrics
                    [sketch_row_space_tranposed,R] = numpy.linalg.qr(sketch.transpose())
                    res['signal_fro_in_sketch_row_space'] = numpy.linalg.norm(numpy.dot(signal_row_space,sketch_row_space_tranposed))
            
                    if not dimension_too_large_for_brute_force:    
                        #diff = ATA - numpy.dot(sketch.transpose(),sketch)
                        diff = A - numpy.dot(A,numpy.dot(sketch.transpose(),sketch))

                        #diff_trace = numpy.trace(diff)
                        #res['diff_trace'] = diff_trace
                        #diff_two_norm = numpy.linalg.norm(diff,2)
                        #res['diff_two_norm'] = diff_two_norm
                        diff_fro_norm = numpy.linalg.norm(diff,'fro')
                        res['diff_fro_norm'] = diff_fro_norm

                        total_res += res['diff_fro_norm']
                        #total_res += res['total_time']

                        value = '%s\n'%json.dumps(res)
                #print 'before write to file'
                        #sys.stdout.write(value)               
                        results_computed += 1

                total_res = total_res/float(5)
                row_to_write.append(total_res)

                
        except:
            traceback.print_exc()
            sys.stderr.write('ERROR IN EXECUTIOM OF MAIN LOOP\n')
    outwriter.writerows(data_to_write)
    outfile.close()
    #book.save("test10.xlsx")
