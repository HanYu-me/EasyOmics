


| | |
| --- | --- |
| Author  | Yu Han ([yuhan](https://www.researchgate.net/profile/Yu-Han-165)) |
| Email   | <yu_han@stu.scu.edu.cn> |

#sad

![Alt text](figures/image.png)
1. Software Installation
2.1 Docker Installation
Docker, an open-source project developed in Go language, packages software and its dependencies into images for quick deployment of the runtime environment. 
For Linux/Mac OS users, the Docker download command:
（1）sudo apt-get install docker-ce docker-ce-cli containerd.io
For Windows system users, follow these steps:
（1）	Refer to Microsoft's documentation to download and set up the WSL2 subsystem: https://docs.microsoft.com/en-us/windows/wsl/install
（2）	Visit Docker's official website to download the Docker Desktop application for Windows: https://www.docker.com/products/docker-desktop
2.2 Software Download
For Linux OS/Mac users, enter the following commands in the terminal：
（1）	docker pull yuhan2000/gwas:1.1.3
（2）	docker images
For Windows users, follow these steps：
（1）	Open the Docker Desktop application.
（2）	Search for and download the image yuhan2000/gwas:1.1.3. 
 
2.3Software Operation and Configuration
For Linux/Mac OS users, use the following commands in the terminal：
（1）	mkdir local_directory 
（2）	docker run -p 6379:6379 -v -p port:3838 local_directory:/srv/shiny-server/Analysis_Result/ yuhan2000/gwas:1.1.3
（3）	Enter ‘localhost:port’ in the browser's address bar to access the software interface.
（4）	docker exec -it gwas 
For Windows users, follow these steps：
（1）	Click to run the downloaded software.
 
（2）	Set the directory and port for output results.
