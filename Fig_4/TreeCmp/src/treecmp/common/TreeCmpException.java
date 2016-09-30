/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package treecmp.common;

/**
 *
 * @author Damian
 */
public class TreeCmpException extends Exception {

  private String errMsg;
//----------------------------------------------
// Default constructor - initializes instance variable to unknown
  public TreeCmpException()
  {
    super();             // call superclass constructor
    errMsg = "unknown";
  }

//-----------------------------------------------
// Constructor receives some kind of message that is saved in an instance variable.
  public TreeCmpException(String errMsg)
  {
    super(errMsg);     // call super class constructor
    this.errMsg = errMsg;  // save message
  }

//------------------------------------------------
// public method, callable by exception catcher. It returns the error message.
  public String getError()
  {
    return errMsg;
  }
}




