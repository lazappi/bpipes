
hello = {
    exec "echo Hello", "hello"
}

world = {
    exec "echo World!", "world"
}

hello_world = {
    exec """
        echo Hello;
        echo World!
    """, "hello_world"
}
