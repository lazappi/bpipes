
hello = {
    exec "echo Hello"
}

world = {
    exec "echo World!"
}

hello_world = {
    exec """
        echo Hello;
        echo World!
    """
}
