

# Objectives
The aim of this first session is to get you familiar with the `fitcourseR` package. This `R` package contains 2 classes and several functions. Let's first have a look at the classes `fitparam` and `fitmodel`. These classes are intended to help you to structure your code and help us to debug it ;)

# Classes
## Fitparam
This class allows you to create a `fitparam` object for each of your model parameters. This object will holds all the properties of your parameter. Some properties are required, like it's name and it's value, whereas others are optional and will be introduced later. These properties are passed as argument to the constructor `fitparam()`, which will return an object of the class `fitparam`: 


```r
R0 <- fitparam(name="R0",value=2)
infectious.period <- fitparam(name="IP",value=5)
class(R0)
```

```
## [1] "fitparam"
```

Now, let's group our `fitparam` objects in a `list`: 

```r
list.fitparams <- list(R0,infectious.period)
```

However, most of the time you will only need to get the values of all your parameters as a named numeric vector `theta`. The function `getParameterValues()` will do that for you:


```r
theta <- getParameterValues(list.fitparams)
print(theta)
```

```
## R0 IP 
##  2  5
```

alternatively, if you want to modify the `value` of one or more `fitparam` objects, just pass a named vector and your list of `fitparam` objects to the function `setParameterValues()`:


```r
new.value <- c(R0=5)
new.list.fitparams <- setParameterValues(list.fitparams,new.value)
getParameterValues(new.list.fitparams)
```

```
## R0 IP 
##  5  5
```

Finally, the `fitparam` constructor performs several checks. For instance, if specify the support of your parameter (optional argument) it will check that the `value` provided is within the support:


```r
R0 <- fitparam(name="R0",value=-1,support=c(0,Inf))
```

```
## Error: 'value' argument is not within the support
```

For more details and examples just type `?fitparam`.

## Fitmodel
In the same way as a `fitparam` object contains all the information related to a parameter, a `fitmodel` contains all you need to fit your model. For instance, it will include the parameters of your model (passed as a list of `fitparam` objects), a function to simulate your model, the data you want to fit your model to, a function to compute the log-likelihood, etc.

You will introduce these elements step by step into your `fitmodel`, which will check that they are consistent and return an error message if not. So please pay attention to the error messages.. if you have any ;)

In the next section, you will get more familiar with this class by building your first `fitmodel` but you can already have a look at the documentation (`?fitmodel`).

# Functions

The package contains 3 types of functions and the fun part will be to re-code most of them:

1. __Functions belonging to a `fitmodel`:__ these are specific to the model and will be passed to `fitmodel` constructor, which will automatically check their arguments, behaviour and returned values. The help of `fitmodel` describes all these functions in more details and illustrates them on a SIR `fitmodel` example.
2. __Functions taking a `fitmodel` as argument:__ these are generic in the sense that they should work for any `fitmodel`. _You should keep this in mind when you'll code them_. They won't be automatically checked so in theory you are free to  specify what they should take as argument and what they should return. That said, we recommend to use the skeleton of the package function for compatibility issues (i.e. if you want to plug your function with another one in the package).
3. __Functions useful and helpful:__ these are short generic functions that are not particularly interesting for you to re-code. Some are like small bricks that you might find useful when you're coding functions of the first two types. Some other are helpful to analyse and visualise your outputs. Please refer to their helps for more details on how to use them.


To re-code the functions of the first two types, we will provide you with their skeleton:


```r
my_greatFunction <- function(arg.1, arg.2){
    ## some guidance about what to code here
}
```

Note that we'll use the `lowerCamelCase` convention for function names and the `dot.separated` convention for argument and variable names. Please try to use the same convention and not to change the names of the arguments for the function of the first category. Of course, since these functions are already available in the package you can use them if you are in trouble or have a quick look how they are coded for comparison with yours. You can access them by removing the `my_` from the skeleton names:


```r
## documentation
?greatFunction
## code
print(greatFunction)
```

# Navigate
Next: [My first fitmodel](first_fitmodel.md)
Previous: [Installation](README_new.md)
