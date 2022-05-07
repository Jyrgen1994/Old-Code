import psycopg2
from psycopg2 import Error
import pandas as pd

class Data:

    # Constructor
    def __init__(self):
        super().__init__()#tror ej denna behövs då vi inte har en "superklass"
        self.connection, self.cursor = self.database_connection()

    #Connecting to database
    def database_connection(self,
                            user = "aissata",
                            host = "127.0.0.1",
                            port = "5432",
                            database = "lms"):
        '''
        ===INITIAL CONNECTION TO DATABASE===
        INPUT: connection infromation for database
        OUTPUT: connection & cursor to database
        '''
        #Try connection to database
        try:
            connection = psycopg2.connect(user = user,
                                          host = host,
                                          port = port,
                                          database = database)

            cursor = connection.cursor()
            # Print PostgreSQL Connection properties
            print ( connection.get_dsn_parameters(),"\n")

            # Print PostgreSQL version
            cursor.execute("SELECT version();")
            record = cursor.fetchone()
            print("You are connected to - ", record,"\n")
            return connection, cursor

        except (Exception, psycopg2.Error) as error :
            print ("Error while connecting to PostgreSQL", error)
            return -1

    #closing database connection.
    def close_connection(self):
        '''
        ===CLOSING CONNECTION===
        '''
        if(self.connection):
            self.cursor.close()
            self.connection.close()
            print("PostgreSQL connection is closed")


    #Insert data to table in database
    def insert_data_table(self, data_frame_processed, table_name):
        '''
        ===INSERTION OF DATA INTO A DATABASE TABLE===
        INPUT: pandas data_frame for processed data
        OUTPUT: None (writes to database)
        '''
        insert_table_query = ' INSERT INTO ' + "\"+" + str(table_name) + "\"" +  ' ('

        #Concatinate attributes to SQL-Query:
        #Iteare over all keys
        for index, attribute in enumerate(data_frame_processed.keys()):
            if index == data_frame_processed.shape[1] -1:
                insert_table_query += attribute
            else:
                insert_table_query += attribute + ', '

        insert_table_query += ')' + ' VALUES \n'

        #Concatinate values to SQL-Query
        #Iterate over rows in data frame
        for index, row in data_frame_processed.iterrows():
            insert_table_query += "("
            #Iterate over keys (columns) in data frame
            for key_index, key in enumerate(data_frame_processed.keys()):
                #END IF FOR EACH ROW
                if key_index == data_frame_processed.shape[1] - 1:

                    #END IF FOR ALL DATA
                    if index == data_frame_processed.shape[0] - 1 :
                        insert_table_query += str(row[key]) + ");"
                    else:
                        insert_table_query += str(row[key]) + ") ,\n"
                else:
                    insert_table_query += str(row[key]) + " , "

        #UNCOMMENT TO SEE how SQL looks, it's human readable
        print(insert_table_query)


        #Try to execute insert
        try:
            self.cursor.execute(insert_table_query)
        except Exception as error:
            print("There was an error iserting data to table : " , error)
        #TO ADD: COMMIT CONDITIONS default is to NOT commit
        if False:
            self.connection.commit()


    def get_data_table(self, table_name):
        '''
        ===FETCHING DATA FROM DATATABLE===
        '''
        select_table_query= 'SELECT * FROM ' + str(table_name) + ';'

        try:
            self.cursor.execute(select_table_query)
        except Exception as error:
            print("There was an error selecting data from table : " , error)
            return -1

        fetched_table = self.cursor.fetchall()
        colnames = [desc[0] for desc in self.cursor.description]

        df = pd.DataFrame(fetched_table, columns=colnames)
        #print(df)
        #print(df.loc[0:2,colnames])
        return df





    # table creation

    # def database_create_table(self, table_name,nr_of_keys,attributes):
    #    '''
    #    NR_OF_KEYS = NR OF PRIMARY KEY
    #    ATTRIBUTES = [(name,data_type),...]
    #    '''

    #    create_table_query= 'CREATE TABLE ' + str(table_name)

    #    '''
    #    (ID INT PRIMARY KEY     NOT NULL,
    #    NAME    VARCHAR (50) NOT NULL);'''

    #    self.cursor.execute(create_table_query)
    #    self.connection.commit()
    #    print("Politician table created!")
