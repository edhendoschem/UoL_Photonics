# Communicating using TCP/IP and Bluetooth

## Table of contents
1. [What is TCP/IP](#item1)
2. [The socket library](#item2)
3. [Client side socket](#item3)
4. [Server side socket](#item4)
5. [Continuous data streaming](#item5)
6. [Bluetooth client and server](#item6)
7. [Further reading](#item7)


## What is TCP/IP? <a name="item1"></a>
TCP/IP is the abbreviation for Transmission Control Protocol/Internet Protocol and is a set of rules specifying how
data should be exchanged. This is the protocol the internet uses and can be leveraged in Python to transmit and receive
useful information and remote control equipment.


## The socket library <a name="item2"></a>
'socket' is a library in Python that allows us to easily create communication pipelines between a client and a server.
Newer versions of Python come with the socket library by default, it allows a number of communication protocols and 
modes. We will focus on the commonly used 'stream' mode for both TCP/IP and BLUETOOTH

In order to create a socket to use TCP/IP protocol we type:
```python
import socket
s = socket.socet(socket.AF_INET, socket.SOCK_STREAM)
```

For a bluetooth socket we use:
```python
import socket
s = socket.socet(socket.AF_BLUETOOTH, socket.SOCK_STREAM, socket.BTPROTO)
```


## Client side socket <a name="item3"></a>
Given an ip address and port, the client side will attempt to connect to a server in listening mode.
We need to find out the ip address of the server computer. To do this in windows we click on start -> run -> cmd
and the type:
```
ipconfig
```
From there we record the value for IPv4 'ddd.ddd.ddd.ddd'  where 'ddd' can be 1, 2 or 3 digits e.g. (10.45.104.2)
We now need a port value, which can be any integer between 1 to 65535. It is recommended to use a higher value such as
5000 to avoid interfering with ports regularly used by other programs.
In order to create the connection with the server we use:
```python
s.connect((ip_address, port))
```
Note that ip_address, port are inside a tuple hence the double '(())', one for the function call and the other for the
tuple.

If the server is listening for connections and we entered all the information correctly then we should be connected
and be able to receive information. It is good practice when setting up the server to send a message to the client
acknowledging the establishment of a connection

<b>Note:</b>
For testing purposes it is useful to have the computer act as both client and server. In order to do this use as ip
address '127.0.0.1' for both client and server.

We can now start sending and receiving data. In order to send data we need to convert it to bytes first, a convenient
way to do this and to improve readability would be to define a function to do it for us:
```python
def packet(obj):
    obj = str(obj)
    obj = bytes(obj, 'utf-8')
    return obj
```

We can now send message like this:
```python
s.send(packet(our_message))
```
If we are interested in sending a string, we can use our previously defined 'packet()' function or we can simply preface
the string with 'b' and it will be automatically converted:
```python
s.send(b'hello!')
```

If we are expecting a response from the server then we need to use the '.recv(n)' method:
```python
data = s.recv(5)
```
The number 'n' used as an argument for '.recv(n)' is the message size in bytes, if the size is smaller than the message then
the message will be truncated and the rest of the message will be sent in the next packet or will be lost if the
client or server closes the connection

Finally, once the connection is no longer needed we close it:
```python
s.close()
```
Once the connection is closed from the client side and the server side, then the port will be
freed and usable again.

The entire python script for the client side is in the file TCP-IP_tut_client file:
```python
#Client side script
import socket #loads the socket library
def packet(obj):
    """
    Takes as argument any object that can be converted to string and converts it to bytes
    """
    obj = str(obj)
    obj = bytes(obj, 'utf-8')
    return obj
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) #Establish the type of socket we want
s.connect(('127.0.0.1', 5003)) #Connect to the desired ip address and port
s.send(b'hello') #Send the message as a bytes object
s.send(packet(102.35))
data = s.recv(5) #Receive a message of 5 bytes in size
print(data) #Print the received message
s.close() #Close the connection
```


## Server side socket <a name="item4"></a>

The server side socket is quite similar to the client side, but instead of using connect, we need to bind and listen 
for connections and then accept the connection with the methods '.listen(n)' and '.accept()':
```python
import socket
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5003))
s.listen(1) #Max number of connections
connection, address = s.accept()
```
The program will remain in the s.listen() until a request to connect has been received. The 's.accept()' line will
return a socket object and a tuple containing ip address and port number. In order to send or receive data, we need
to interact with the socket object, therefore we store it in the variable 'connection'.

To send or receive data is as straightforward as before:
```python
connection.send(b"Connection established")
received_data = connection.recv(32) #received data is stored in the variable 'data'
```

Once we have finished receiving data we can close the connection (if the client hasn't) and then close the socket by
using:
```python
connection.close()
s.close()
```

The complete file TCP-IP_tut_server should look like this:
```python
import socket
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5003))
s.listen(1)
connection, address = s.accept()
received_data = connection.recv(32)
connection.send(b"Connection established")
print(received_data) #
connection.close()
s.close()
```


## Continuous data streaming <a name="item5"></a>
In our previous client/server examples, the client and server would each send and receive exactly one packet of data
which in most cases isn't very useful. In order to have the server continuously be able to send and receive data, we need to
place our send/receive commands inside endless loops. For the server side we have:

```python
import socket #loads the socket library
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5003))
s.listen(1)
connection, address = s.accept()
connection.send(b"Connection established")
while True: #Endless loop
    received_data = connection.recv(32) #Receives packets of 32 bytes in size
    print(received_data) #Prints the data received (optional)
    if not received_data: #If the connection was closed or an empty packet sent, breaks the endless loop
        break;
    connection.send(b"data received") #Sends a message to the client to prompt the next packet
connection.close()
s.close()
```

The previous code allow us to receive endless amount of data and continuously store it in the 'received_data' variable
to be used however we see fit, before it is overwritten in the next cycle. Since we are not creating new variables
in each cycle, then the amount of RAM memory consumed by the script should remain constant. At the end of the cycle
a message is sent to the client to unpause its process and prompt the next data packet.

If the client closes the connection then an empty packet is sent by dedaulf which will trigger the break clause, 
ending the loop.

We can use the data to perform any action inside the loop including calling other functions, keeping in mind that the 
data received is in bytes and will need to be converted to a string or number depending on what we need. Here is an 
example of writing the data being sent by the client as .csv (comma separated values) to a file:
```python
import socket #loads the socket library
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5003))
s.listen(1)
connection, address = s.accept()
connection.send(b"Connection established")
file_to_write = open("test_output.csv", "a+") #Opens a file and if it doesn't exist creates it, in append mode
while True:
    received_data = connection.recv(32)
    print(received_data)
    if not received_data:
        break;
    #Converts the bytes data to a string, removes the "b" and "'", adds a newline character and stores the result 
    #at the end of the file
    file_to_write.write(str(received_data).lstrip('b').strip("'")+'\n') 
    connection.send(b"data processed")
file_to_write.close() #Closes the file
connection.close()
s.close()
```

If we are to receive large volumes of data, we might want to avoid the line:
```python
print(received_data)
```
as this will continuously print what we receive in a new line until the computer runs out of memory. Another way to
solve this issue would be to use the optional argument for the 'print' function 'end = "\r"' which will
overwrite the same line, therefore preventing an endless increase of memory usage. Our final server file looks like 
this:
```python
#Improved Server side script
import socket #loads the socket library
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5004))
s.listen(1)
connection, address = s.accept()
connection.send(b"Connection established")
file_to_write = open("some_file.csv", "a+")
while True:
    received_data = connection.recv(128)
    print(received_data, end = '\r') #Will print the data received and continuously overwrite this line
    if not received_data:
        break;
    #Converts the bytes data to a string, removes the "b" and "'", adds a newline character and stores the result 
    #at the end of the file
    file_to_write.write(str(received_data).lstrip('b').strip("'")+'\n')
    connection.send(b"data processed")
file_to_write.close()
connection.close()
s.close()
```

Now let's say we have a datalog with some measurements being continuously updated, we want to send that data as it
is being written to our 'server' computer for processing and/or storage. We could try opening the file and then
setting an endless loop to send the data, however this approach has two problems:

First, the file opened will be a copy of the file being updated by whatever measurement program we're interested in.
This means, it will contain all the lines up to the point the 'open' command was executed. Any new lines written after
that will not be included.

Second, even if we create an endless loop to constantly reopen the file to obtain the new data being written, we will
eventually run out of memory as the file gets progressively larger, not to mention that opening the file will become
progressively slower.

To deal with these problems we need to use a generator, a function that returns an iterator. The main characteristic of
iterators is that they are only 'aware' of the current element and store no previous elements, which makes their memory
usage small and constant. We can explicitly call the next element on an iterator using the 'next(iterator_name_here)' 
function, or implicitly by using a 'for' loop. In order to define a function that returns an iterator we need to
end it in 'yield' instead of 'return'.

To create a function that returns an iterator of lines in a file we use:
```python
def line_get(file_handle):
    """
    Returns a generator that endlessly reads a file while occupying constant memory. 
    requires file_handle(file object)
    """
    file_handle.seek(0,1) #Change this from (0,1) to (0,2) if you want to start reading from the last line
    while True:
        line = file_handle.readline()
        if not line:
            time.sleep(10) #Waits 10 seconds before attempting to read the next line
            continue
        yield line
```

This function will return an iterator that yields lines on the file we used as argument, each time we call
'next()' explicitly or implicitly. If it cannot find more lines, waits 10 seconds until attempting to read the
next line. This function will read the log file from the beginning, but if we are only interested in the most recent
value, then we can set:
```python
file_handle.seek(0,2)
```
This makes the function start at the end of the current file, which means only any new lines written after that will
be read.

With that function we can now define our client to read and transmit data from a file as it is being written. To
simulate a file being written we have provided data_logging.ipynb which continuously appends data to a the file
'testfile.csv' in 5 second intervals:
```python
import time
import numpy as np

file_handle = open("testfile.csv", 'a+', 1) #Opens in append mode to add to the end of the file, 
					    #flushes (copies to file) after every line
try:
    while True: #Needs to be stopped manually by clicking or selecting interrupt kernel in jupyter notebook
        rand1 = str(np.random.randint(0,100000)) #Random data
        rand2 = str(np.random.randint(0,100000)) #Random data
        rand3 = str(np.random.randint(0,100000)) #Random data
        file_handle.write(rand1+','+rand2+','+rand3+'\n') #write the random data
        time.sleep(5) #wait 5 seconds
except: #will catch the manual stop and close the file
    file_handle.close()
```

The client can be defined as:
```python
#Improved client side script, endlessly read a log file
import socket #loads the socket library
import time

def packet(obj):
    """
    Takes as argument any object that can be converted to string and converts it to bytes
    """
    obj = str(obj)
    obj = bytes(obj, 'utf-8')
    return obj

def line_get(file_handle):
    """
    Returns a generator that endlessly reads a file while occupying constant memory. 
    requires file_handle(file object)
    """
    file_handle.seek(0,1) #Change this to (0,2) if you want to start reading from the last line
    while True:
        line = file_handle.readline()
        if not line:
            time.sleep(10) #Waits 10 seconds before attempting to read the next line
            continue
        yield line
    

file_to_read = open('testfile.csv', 'r')
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) #Establish the type of socket we want
s.connect(('127.0.0.1', 5004)) #Connect to the desired ip address and port
print(s.recv(128)) #Prints connection established message from server
try:
    for lines in line_get(file_to_read):#Needs to be manually interrupted using the stop button in jupyter notebook
        s.send(packet(lines)) #Send the message as a bytes object
        data = s.recv(128) #Receive a message, cycle stops here until the message is received
        print(data, end = '\r') #Print the received message, will continuously overwrite this line
except: #Catches the keyboard interrupt and closes the connection
    s.close() #Close the connection
```

This code is similar to our original client code, however the 'try: except:' clause has been added to catch the 
KeyboardInterrupt error we get when manually stopping our endless iterator, making sure it closes the connection.
Successfully closing the connection is important to signal our server to stop listening and close the connection.

The 'for lines in line_get(file_to_read):' line, endlessly calls on our generator function to provide the lines read
from the file 'testfile.csv'. An explicit way to write this would be:
```python
gen = line_get(file_read) #creates the iterator
while True:
    line = next(gen)
```


## Bluetooth client and server <a name="item6"></a>
In order to create a client and server in bluetooth we can reuse the same code with the only change being:
```python
s = socket.socket(socket.AF_BLUETOOTH, socket.SOCK_STREAM, socket.BTPROTO_RFCOMM)
```
And instead of using ip addresses, we need to know the MAC address of the server in the form 'nn:nn:nn:nn:nn:nn', 
where n can be either a number or a letter e.g 78:F8:2C:F0:2B:EC, and the port can be set to any number between 1 to 30

<b>Note</b>: Bluetooth socket works fine in linux, but on Windows it doesn't seem to work. The module socket.AF_BLUETOOTH is
missing. An alternative could be the pybluez library, but it isn't currently available for Python version 3.6 
(14/09/2017)


## Further reading <a name="item7"></a>
Sockets:  
https://docs.python.org/3/howto/sockets.html  
http://blog.kevindoran.co/bluetooth-programming-with-python-3/

Iterators and generators:  
https://docs.python.org/3/howto/functional.html
