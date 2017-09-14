# Communicating using TCP/IP and Bluetooth

## Table of contents
1. [What is TCP/IP](#item1)
2. [The socket library](#item2)
3. [Client side socket](#item3)
4. [Server side socket](#item4)
5. [Continuous data streaming](#item5)
6. [Bluetooth client and server](#item6)


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
From there we record the value for IPv4 ddd.ddd.ddd.ddd  where ddd can be 1, 2 or 3 digits e.g. (10.45.104.2)
We now need a port value, which can be any integer between 1 to 65535. It is recommended to use a higher value such as
5000 to avoid interfering with ports regularly used by other programs.
In order to create the connection with the server we use:
```python
s.connect((ip_address, port))
```
Note that ip_address, port are inside a tuple hence the double (()), one for the function call and the other for the
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

If we are expecting a response from the server then we need to use the '.recv()' method:
```python
data = s.recv(5)
```
The number used as an argument for '.recv()' is the message size in bytes, if the size is smaller than the message then
the message will be truncated and the rest of the message will be sent in the next packet or will be lost if the
client or server closes the connection

Finally, once the connection is no longer needed we close it:
```python
s.close()
```
Once the connection is closed from the client side and if successfully from the server side, then the port will be
freed and usable again

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
for connections and then accept the connection with the methods '.listen()' and '.accept()':
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
which in most cases isn't very useful. In order to have the server continuosly be able and receive data, we need to
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
in each cycle then the amount of RAM memory consumed by the script should remain constant. At the end of the cycle
a message is sent to the client to unpause its process and prompt the next data packet.

If the client closes the connection then an empty packet will be sent which will trigger the break clause, ending the
loop

We can use the data to perform any action inside the loop including calling other functions, keep in mind that the data
received is in bytes and will need to be converted to a string or number depending on what we need. Here is an example
of writing the data being sent by the client to a file:
```python
import socket #loads the socket library
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 5003))
s.listen(1)
connection, address = s.accept()
connection.send(b"Connection established")
file = open("some_file.csv", "a+") #Opens a file and if it doesn't exist creates it, in append mode
while True:
    received_data = connection.recv(32)
    print(received_data)
    if not received_data:
        break;
    #Converts the bytes data to a string, removes the "b" and "'", adds a newline character and stores the result 
    #at the end of the file
    file.write(str(received_data).lstrip('b').strip("'")+'\n') 
    connection.send(b"data processed")
file.close() #Closes the file
connection.close()
s.close()
```



